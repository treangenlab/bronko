use crate::cli::*;
use crate::lcb::*;
use crate::util::*;
use crate::build::*;

use anyhow::{Result};

use num_cpus;
use log::*;
use crate::consts::*;
use rustc_hash::{FxHashMap};
use dashmap::DashMap;
use bincode::{config};

use std::fs::File;
use std::vec;
use std::fs;

use rayon::prelude::*;
use rayon::join;

use std::path::{Path};
use std::io::{BufReader, BufRead, BufWriter, Write};
use std::process::{Command, Stdio};


fn check_args(args: &CallArgs) {
    let output_level;
    if args.verbose {
        output_level = log::LevelFilter::Trace;
    } else if args.debug {
        output_level = log::LevelFilter::Debug;
    } else {
        output_level = log::LevelFilter::Info;
    }

    simple_logger::SimpleLogger::new()
        .with_level(output_level)
        .init()
        .unwrap();

    //Check kmer size to make sure it is odd and greater than 3
    if args.kmer % 2 != 1 || args.kmer > MAX_KMER_SIZE || args.kmer < MIN_KMER_SIZE {
        error!("Invalid kmer size, must be odd and between [{}-{}]", MIN_KMER_SIZE, MAX_KMER_SIZE);
        std::process::exit(1)
    }

    //check to see if inputs are valid fastq and fasta files
    for fastq_file in &args.reads {
        if !check_fastq(&fastq_file) {
            error!("{} does not appear to be a fastq file (must be .fq(.gz)/.fastq(.gz)/.fnq(.gz))", &fastq_file);
            std::process::exit(1)
        }
    }

    if args.genomes.is_some() && args.db.is_some() {
        error!("Please provide either a db or the genomes you would like to index, not both.");
        std::process::exit(1);
    } else if args.genomes.is_none() && args.db.is_none() {
        error!("Please provide either a db or the genomes you would like to index.");
        std::process::exit(1);
    }

    if let Some(genomes) =  &args.genomes {
        for fasta_file in genomes {
            if !check_fasta(&fasta_file){
                error!("{} does not appear to be a fasta file (must be .fa(.gz)/.fasta(.gz)/.fna(.gz))", &fasta_file);
                std::process::exit(1)
            }
        }
    }

    let available_threads = num_cpus::get();
    if args.threads <= 0 {
        error!("Number of threads must be greater than 0");
        std::process::exit(1);
    } else if args.threads > available_threads {
        error!("You requested {} threads but only have {} available on your system", args.threads, available_threads);
        std::process::exit(1);
    }

    if args.min_af < 0.01 {
        warn!("Minimum allele frequency set below 0.01, more false positive variants will be returned. We suggest setting this to a more realistic threshold (0.01-0.05)");
    } else if args.min_af > 1.0 {
        error!("Minimum allele frequency set above 1, please set between 0-1 (recommended between 0.01-0.05)");
        std::process::exit(1);
    } else if args.min_af >= 0.5 {
        warn!("Minimum allele frequency set equal to or greater than 0.5, no minor variants will be returned");
    }

    if args.n_per_strand <= 0 {
        warn!("Number of kmers per strand set to 0, this is equivalent to no strand filtering")
    } else if args.n_per_strand >= args.kmer {
        error!("Number of kmers per strand set >= k, please set lower value (recommended 2-4, default 2)");
        std::process::exit(1);
    } else if args.n_per_strand >= 5 {
        warn!("Number of kmers per strand set very high, only strongly supported variants will be returned")
    }

    if args.first_pairs.len() != args.second_pairs.len() {
        error!("Number of paired end sequences do not match, exiting.");
        std::process::exit(1);
    }

}

#[derive(Debug)]
pub struct OutputInfo {
    filename: String,
    selected_genome: String,
    num_major_variants: usize,
    num_minor_variants: usize,
    breadth_coverage: f64,
    depth_coverage: f64,
    num_perfect_kmers: usize,
    num_variant_kmers: usize,
    num_unmapped_kmers: usize
}

pub fn call(args: CallArgs) {

    //First check to make sure none of the arguments are invalid
    check_args(&args);
    trace!("k={}, threads={}", args.kmer, args.threads);

    // create output directory
    let out_path = Path::new(&args.output);
    fs::create_dir_all(out_path).unwrap_or_else(|e| {
        error!("{} | Unable to create outputs in output directory 2", e);
        std::process::exit(1);
    });



    let (ref_index, viral_metadata);

    
    //build the indexes
    if let Some(genomes) = &args.genomes {
        info!("Creating bronko index from provided reference genomes");
        let (index, meta): (FxHashMap<u64, Vec<BucketInfo>>, ViralMetadata) = build_indexes(args.kmer, &genomes).unwrap_or_else(|e| {
            error!("{} | Reference failed to build", e);
            std::process::exit(1)
        });
        ref_index = index;
        viral_metadata = meta;
        log_memory_usage(true, "Fasta files indexed successfully. Starting counting kmers ");
    } else if let Some(db) = &args.db {
        info!("Reading in provided bronko index");
        let file = File::open(&db).unwrap_or_else(|e| {
            error!("Failed to open file '{}': {}", &db, e);
            std::process::exit(1);
        });

        let mut reader = BufReader::new(file);
        let config = config::standard();
        let db: BronkoIndex = bincode::decode_from_std_read(&mut reader, config).unwrap_or_else(|e| {
            error!("Failed to read Bronko Index from '{}': {}", db, e);
            std::process::exit(1);
        });

        let db_k = db.k;
        if db_k != args.kmer {
            error!("Database k is not the same as provided, please set -k to {} or build a new index", {db_k});
            std::process::exit(1);
        }
        
        ref_index = db.global_index;
        viral_metadata = db.metadata;
    } else {
        error!("Unable to build/read index, exiting");
        std::process::exit(1);
    }

    // storing the variant information
    let total_samples = args.reads.len() + args.first_pairs.len();
    let mut variant_info: Vec<(String, Vec<VCFRecord>)> = Vec::with_capacity(total_samples); // (sample_name, variant records)
    let mut output_info: Vec<OutputInfo> = Vec::with_capacity(total_samples); //storing the outputs for tsv overview

    // PROCESS SINGLE END READS
    if args.reads.len() > 0 {
        for se_read in args.reads.iter() {
            info!("Processing {}", se_read);

            //Get kmer counts
            let (kmers, total_reads, total_kmers, unique_kmers, unique_counted_kmer ) = get_kmers(&se_read, &args.threads, &args);
            info!("{} reads counted from {}", total_reads, se_read);
            info!("{} unique kmers above {} count, {} total unique kmers, {} total kmers (~{} basepairs)", unique_counted_kmer, args.min_kmers, unique_kmers, total_kmers, total_kmers*args.kmer);
            log_memory_usage(true, "Finished counting kmers");

            //initialize output storage and then map the kmers using the index
            info!("Initializing mapping arrays");
            let output_maps = initialize_output_maps(&viral_metadata);
            info!("Mapping kmers to all genomes");
            let mapping_data = map_kmers(&kmers, &ref_index, &viral_metadata, &args.threads, &args, &output_maps);

            //select the best genome from the mapping data
            info!("Selecting the most representative genome");
            let best_genome_index = pick_best_genome(&mapping_data, &viral_metadata).unwrap_or_else(|| {
                error!("Unable to pick a best genome");
                std::process::exit(1);
            });

            //print out the best selected genome and mapping statistics
            let (n_perfect_mapped, n_variant_mapped, n_unique_perfect) = mapping_data.get(&best_genome_index).unwrap_or_else(|| {
                error!("Error getting mapping statistics for best genome");
                std::process::exit(1);
            });
            let best_genome_filename = &viral_metadata.files[best_genome_index as usize].name;
            info!("Selected a representative genome: {}", best_genome_filename);
            let n_unmapped_kmers = unique_counted_kmer-n_perfect_mapped-n_variant_mapped;
            let message = format!("Mapped {}/{} kmers perfectly ({} unique among refs), {}/{} had a variant, {} unmapped", n_perfect_mapped, unique_counted_kmer, n_unique_perfect, n_variant_mapped, unique_counted_kmer, n_unmapped_kmers); 
            log_memory_usage(true, &message);

            if ((n_variant_mapped + n_perfect_mapped) as f64 / unique_counted_kmer as f64) < 0.2 {
                warn!("Percent of kmers found is very low for this reference, suggesting lack of a representative reference, a bad sequencing run, contamination in sample, or some other issue")
            }

            
            // call variants and print them out to vcf
            let (output, output_rev, output_counts, output_rev_counts)= output_maps.get(&best_genome_index).unwrap_or_else(|| {
                error!("Failed to find mapping data for selected genome");
                std::process::exit(1);
            });

            let (variants, num_major_variants, num_minor_variants, breadth_cov, depth_cov) = call_variants(
                &args,
                output,
                output_counts,
                output_rev,
                output_rev_counts,
                args.min_af,
                !args.no_end_filter,
                !args.no_strand_filter,
                args.n_per_strand,
                &best_genome_filename
            );
            log_memory_usage(true, "Called variants successfully");

            //print outputs
            if args.output_pileup {
                print_pileup(&se_read, &args, &output, &output_rev, &viral_metadata, &best_genome_index);
            }
            print_output(&se_read, &args, &variants, &viral_metadata, &best_genome_index);

            output_info.push(OutputInfo {
                filename: se_read.to_string(),
                selected_genome: best_genome_filename.to_string(),
                num_major_variants: num_major_variants,
                num_minor_variants: num_minor_variants,
                breadth_coverage: breadth_cov,
                depth_coverage: depth_cov,
                num_perfect_kmers: *n_perfect_mapped,
                num_variant_kmers: *n_variant_mapped,
                num_unmapped_kmers: n_unmapped_kmers,
            });
            variant_info.push((se_read.to_string(), variants));

            if !args.keep_kmer_counts {
                cleanup_kmc_files(&args.output);
            }
        }
    }

    // PROCESS PAIRED END READS
    if args.first_pairs.len() > 0 && args.second_pairs.len() > 0 {
        for (r1, r2) in args.first_pairs.iter().zip(args.second_pairs.iter()){
            info!("Processing paired reads {}, {}", r1, r2);

            let half_threads = &args.threads / 2;
            let ((kmers1, total_reads_r1, total_kmers_r1, unique_kmers_r1, unique_counted_kmer_r1),
                 (kmers2, total_reads_r2, total_kmers_r2, unique_kmers_r2, unique_counted_kmer_r2)) =
                join(
                    || get_kmers(r1, &half_threads, &args),
                    || get_kmers(r2, &half_threads, &args),
                );
            info!("{} reads counted from {}", total_reads_r1 + total_reads_r2, r1);
            info!("{} unique kmers above {} count, {} total unique kmers, {} total kmers (~{} basepairs)", unique_counted_kmer_r1 + unique_counted_kmer_r2, args.min_kmers, unique_kmers_r1 + unique_kmers_r2, total_kmers_r1 + total_kmers_r2, (total_kmers_r1 + total_kmers_r2) * args.kmer);
            log_memory_usage(true, "Finished counting kmers");

            //initialize output storage and then map the kmers using the index
            info!("Initializing mapping arrays");
            let output_maps = initialize_output_maps(&viral_metadata);
            info!("Mapping kmers to all genomes");
            let mapping_data_r1 = map_kmers(&kmers1, &ref_index, &viral_metadata, &half_threads, &args, &output_maps);
            let mapping_data_r2 = map_kmers(&kmers2, &ref_index, &viral_metadata, &half_threads, &args, &output_maps);

            info!("Selecting the most representative genome");
            let best_genome_index = pick_best_genome_paired(&mapping_data_r1, &mapping_data_r2, &viral_metadata).unwrap_or_else(|| {
                error!("Unable to pick a best genome");
                std::process::exit(1);
            });

            let (n_perfect_mapped_r1, n_variant_mapped_r1, n_unique_perfect_r1) = mapping_data_r1.get(&best_genome_index).unwrap_or_else(|| {
                error!("Error getting mapping statistics for best genome");
                std::process::exit(1);
            });
            let (n_perfect_mapped_r2, n_variant_mapped_r2, n_unique_perfect_r2) = mapping_data_r2.get(&best_genome_index).unwrap_or_else(|| {
                error!("Error getting mapping statistics for best genome");
                std::process::exit(1);
            });
            
            let best_genome_filename = &viral_metadata.files[best_genome_index as usize].name;
            info!("Selected a representative genome: {}", best_genome_filename);
            let n_unmapped_kmers = (unique_counted_kmer_r1 + unique_counted_kmer_r2)-(n_perfect_mapped_r1+n_perfect_mapped_r2)-(n_variant_mapped_r1 + n_variant_mapped_r2);
            let message = format!("Mapped {}/{} kmers perfectly ({} unique among references) , {}/{} had a variant, {} unmapped", n_perfect_mapped_r1 + n_perfect_mapped_r2, unique_counted_kmer_r1 + unique_counted_kmer_r2, n_unique_perfect_r1 + n_unique_perfect_r2, n_variant_mapped_r1+n_variant_mapped_r2, unique_counted_kmer_r1 + unique_counted_kmer_r2, n_unmapped_kmers); 
            log_memory_usage(true, &message);
            if ((n_variant_mapped_r1 + n_variant_mapped_r2 + n_perfect_mapped_r1 + n_perfect_mapped_r2) as f64 / ((unique_counted_kmer_r1 as f64) + (unique_counted_kmer_r2 as f64))) < 0.2 {
                warn!("Percent of kmers found is very low, suggesting a bad reference, a bad sequencing run, contamination in sample, or some other issue")
            }

            // call variants and print them out
            let (output, output_rev, output_counts, output_rev_counts)= output_maps.get(&best_genome_index).unwrap_or_else(|| {
                error!("Failed to find mapping data for selected genome");
                std::process::exit(1);
            });

            let (variants, num_major_variants, num_minor_variants, breadth_cov, depth_cov) = call_variants(
                &args,
                output,
                output_counts,
                output_rev,
                output_rev_counts,
                args.min_af,
                !args.no_end_filter,
                !args.no_strand_filter,
                args.n_per_strand,
                &best_genome_filename,
            );
            log_memory_usage(true, "Called variants successfully");

            //print outputs
            if args.output_pileup {
                print_pileup(&r1, &args, &output, &output_rev, &viral_metadata, &best_genome_index);
            }
            print_output(&r1, &args, &variants, &viral_metadata, &best_genome_index);

            output_info.push(OutputInfo {
                filename: r1.to_string(),
                selected_genome: best_genome_filename.to_string(),
                num_major_variants: num_major_variants,
                num_minor_variants: num_minor_variants,
                breadth_coverage: breadth_cov,
                depth_coverage: depth_cov,
                num_perfect_kmers: *n_perfect_mapped_r1 + *n_perfect_mapped_r2,
                num_variant_kmers: *n_variant_mapped_r1 + *n_variant_mapped_r2,
                num_unmapped_kmers: n_unmapped_kmers,
            });
            variant_info.push((r1.to_string(), variants));
            
            
            if !args.keep_kmer_counts {
                cleanup_kmc_files(&args.output);
            }
        }
    }

    info!("Printing overview");
    print_output_info(&args, &output_info);
    info!("All samples processed successfully");

    //Printing out alignments
    if args.output_alignment {
        info!("Building alignment(s)");
        build_alignments_for_genomes(&output_info, &variant_info, &viral_metadata, &args);
    }

    info!("");
    info!("bronko complete!");

}

fn cleanup_kmc_files(output_dir: &str) {
    if let Ok(entries) = fs::read_dir(output_dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
                if name.ends_with("kmc_pre")
                    || name.ends_with("kmc_suf")
                    || name.ends_with("_counts.txt")
                {
                    if let Err(e) = fs::remove_file(&path) {
                        eprintln!("Failed to delete {:?}: {}", path, e);
                    }
                }
            }
        }
    }
}

pub fn pick_best_genome(mapping_data: &FxHashMap<u16, (usize, usize, usize)>, viral_metadata: &ViralMetadata) -> Option<u16> {
    //input data has key of usize that is the index in the Viral Metadata, and value of tuple (usize, usize) with the number of perfectly mapped kmers and variant mapped kmers
    // Value tuple: (perfect_kmers, variant_kmers)
    let mut best_genome: Option<u16> = None;
    let mut best_score: f64 = 0.0;

    for (file_index, (perfect, variant, unique_perfect)) in mapping_data.iter() {
        let genome_len: usize = viral_metadata.files[*file_index as usize]
            .sequences
            .iter()
            .map(|s| s.len)
            .sum();

        let score = (*perfect as f64) / (genome_len as f64) / 2.0;

        let genome_name: String = viral_metadata.files[*file_index as usize].name.clone();
        trace!(
            "Genome {}: perfect={}, variant={}, unique={}, len={}, score={:.4}",
            genome_name, perfect, variant, unique_perfect, genome_len, score
        );

        if score > best_score {
            best_score = score;
            best_genome = Some(*file_index);
        }
    }

    best_genome
}

pub fn pick_best_genome_paired(
    mapping_data1: &FxHashMap<u16, (usize, usize, usize)>,
    mapping_data2: &FxHashMap<u16, (usize, usize, usize)>,
    viral_metadata: &ViralMetadata,
) -> Option<u16> {
    let mut combined: FxHashMap<u16, (usize, usize, usize)> = FxHashMap::default();

    // Combine both mapping datasets
    for (&k, (p, v, u)) in mapping_data1.iter() {
        combined.insert(k, (*p, *v, *u));
    }

    // Add everything from mapping_data2
    for (&k, (p, v, u)) in mapping_data2.iter() {
        combined
            .entry(k)
            .and_modify(|(cp, cv, cu)| {
                *cp += p;
                *cv += v;
                *cu += u;
            })
            .or_insert((*p, *v, *u));
    }

    //get best genome just like for single end reads
    let mut best_genome: Option<u16> = None;
    let mut best_score: f64 = 0.0;

    for (file_index, (perfect, variant, unique_perfect)) in combined.iter() {
        let genome_len: usize = viral_metadata.files[*file_index as usize]
            .sequences
            .iter()
            .map(|s| s.len)
            .sum();

        let score = (*perfect as f64) / (genome_len as f64) / 2.0; //divided by 2.0 bc of fwd and reverse

        let genome_name: String = viral_metadata.files[*file_index as usize].name.clone();
        trace!(
            "Genome {} (paired): perfect={}, variant={}, unique_perfect={}, len={}, score={:.4}",
            genome_name, perfect, variant, unique_perfect, genome_len, score
        );

        if score > best_score {
            best_score = score;
            best_genome = Some(*file_index);
        }
    }

    best_genome
}

pub fn build_alignments_for_genomes(
    output_info: &Vec<OutputInfo>,
    variant_info: &Vec<(String, Vec<VCFRecord>)>,
    viral_metadata: &ViralMetadata,
    args: &CallArgs,
) {
    // Build fast lookup for variants
    let mut variant_map: FxHashMap<&str, &Vec<VCFRecord>> = FxHashMap::default();
    for (fname, vars) in variant_info {
        variant_map.insert(fname.as_str(), vars);
    }

    // Group by selected_genome (the best genomes for each samples)
    let mut genome_map: FxHashMap<String, Vec<(String, Vec<VCFRecord>)>> = FxHashMap::default();

    // loop through the output_info and append from the variant_map
    for oi in output_info {

        if oi.breadth_coverage < 0.90 { // don't add genomes with low coverage to the alignment
            info!("Skipping {} (breadth of coverage = {})", oi.filename, oi.breadth_coverage);
            continue;
        }

        if let Some(vars) = variant_map.get(oi.filename.as_str()) {
            genome_map
                .entry(oi.selected_genome.clone())
                .or_default()
                .push((oi.filename.clone(), (*vars).clone()));
        } else {
            warn!("No variant info found for {}", oi.filename);
        }
    }

    // For each genome with ≥3 samples, build alignment
    for (genome_name, samples) in genome_map {
        if samples.len() < 3 {
            info!("Skipping {} (only {} samples)", genome_name, samples.len());
            continue;
        }

        // Find genome metadata in viral_metadata
        if let Some(file_meta) = viral_metadata.files.iter().find(|f| f.name == genome_name) {
            info!(
                "Building alignment for genome {} with {} samples",
                genome_name,
                samples.len()
            );
            build_alignment_fasta(samples, args, file_meta);
        } else {
            warn!("Genome {} not found in metadata, skipping", genome_name);
        }
    }
}



pub fn build_alignment_fasta(
    sample_variants: Vec<(String, Vec<VCFRecord>)>,
    args: &CallArgs,
    file_meta: &FileMeta,
) {
        
    //first collect all sequence/position pairs with a variant and their reference base, as well as a local version for each sample
    let mut all_positions: FxHashMap<(String, usize), u8> = FxHashMap::default();
    let mut sample_positions: FxHashMap<String, FxHashMap<(String, usize), u8>> = FxHashMap::default();

    for (sample, vcf_records) in &sample_variants {
        sample_positions.insert(sample.clone(), FxHashMap::default());
        for variant in vcf_records {
            if variant.af >= 0.5 { //store the major variants in the full set and the sample set
                all_positions.insert((variant.seq.clone(), variant.pos), variant.ref_base);
                if let Some(sample_map) = sample_positions.get_mut(&sample.clone()) {
                    sample_map.insert((variant.seq.clone(), variant.pos), variant.alt_base);
                }
            }
        }
    }

    // then sort all global positions such that they are ordered by seq and position
    let mut positions: Vec<(String, usize)> = all_positions.keys().cloned().collect();
    positions.sort_unstable();

    // now just loop through each position for each sample, and each time output a string if that variant is present (or the reference if not)
    let fasta_out = File::create(format!("{}/{}.mfa", args.output, file_meta.name)).unwrap_or_else(|e| {
        error!("{} | Failed to create mfa alignment file", e);
        std::process::exit(1);
    });
    let mut writer = BufWriter::new(fasta_out);

    // Output the reference sequence first
    let mut ref_seq = String::with_capacity(positions.len());
    for pos_key in &positions {
        let ref_base = all_positions
            .get(pos_key)
            .map(|&b| nucleotide_bits_to_char(b as u64))
            .unwrap_or('N');
        ref_seq.push(ref_base);
    }
    writeln!(writer, ">{}", file_meta.name).unwrap();
    writeln!(writer, "{}", ref_seq).unwrap();

    for (sample_name, sample_map) in &sample_positions {
        let mut seq = String::with_capacity(positions.len());

        for pos_key in &positions {
            // If the sample has a variant at this position, use the alt_base
            if let Some(&alt_base) = sample_map.get(pos_key) {
                seq.push(nucleotide_bits_to_char(alt_base as u64)); // u8 → char
            } else {
                // Otherwise fall back to reference base from all_positions
                let ref_base = all_positions
                    .get(pos_key)
                    .map(|&b| nucleotide_bits_to_char(b as u64))
                    .unwrap_or('N');
                seq.push(ref_base);
            }
        }

        let sample_out = clean_sample_id(&sample_name);

        writeln!(writer, ">{}", sample_out).unwrap();
        writeln!(writer, "{}", seq).unwrap();
    }

}

pub fn get_kmers(reads_file: &String, threads: &usize, args:&CallArgs) -> (Vec<(String, u64)>, usize, usize, usize, usize){
    //count kmers using kmc3
    let kmc_result = count_kmers_kmc(&reads_file, &threads, &args);
    let (total_reads, total_kmers, unique_kmers, unique_counted_kmer) = match kmc_result {
        Ok(values) => values,
        Err(err_msg) => {
            error!("{}", err_msg);
            std::process::exit(1);
        }
    };


    let file_stem = clean_sample_id(&reads_file);
    let kmers = load_kmers(&format!("{}/{}_counts.txt", args.output, file_stem));

    (kmers, total_reads, total_kmers, unique_kmers, unique_counted_kmer)
}

pub fn print_pileup(
    read_output: &String,
    args: &CallArgs,
    output: &DashMap<String, OutputData>,
    output_rev: &DashMap<String, OutputData>,
    viral_metadata: &ViralMetadata,
    best_genome_index: &u16, 
){
    info!("Writing output to pileup");

    let file_stem = clean_sample_id(read_output);

    let tsv_pileup = File::create(format!("{}/{}.tsv", args.output, file_stem)).unwrap_or_else(|e| {
        error!("{} | Failed to create tsv pileup file", e);
        std::process::exit(1);
    });
    let mut writer = BufWriter::new(tsv_pileup);

    writeln!(writer, "reference\tindex\tref\tA\tC\tG\tT\ta\tc\tg\tt").unwrap();


    let file_meta = &viral_metadata.files[*best_genome_index as usize];

    for seq_entry in &file_meta.sequences {
        let seq = &seq_entry.name;

        let fwd = output.get(seq).expect("Could not match seq to fwd counts");
        let rev = output_rev.get(seq).expect("Could not match seq to rev counts");

        for (i, ref_base) in fwd.ref_bases.iter().enumerate() {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                seq,
                i + 1,
                *ref_base as char,
                fwd.counts[i][0],
                fwd.counts[i][1],
                fwd.counts[i][2],
                fwd.counts[i][3],
                rev.counts[i][0],
                rev.counts[i][1],
                rev.counts[i][2],
                rev.counts[i][3]
            ).unwrap();
        }
    }
}


pub fn print_output_info(args: &CallArgs, output_info: &Vec<OutputInfo>) {
    let file_stem = "bronko_overview";
    let output = &args.output;
    let path = format!("{}/{}.tsv", output, file_stem);

    let tsv_file = File::create(&path).unwrap_or_else(|e| {
        error!("{} | Failed to create tsv file", e);
        std::process::exit(1);
    });

    let mut writer = BufWriter::new(tsv_file);

    // header
    writeln!(
        writer,
        "filename\tselected_genome\tnum_major_variants\tnum_minor_variants\tbreadth_coverage\tdepth_coverage\tnum_perfect_kmers\tnum_variant_kmers\tnum_unmapped_kmers"
    ).unwrap();

    // rows
    for info in output_info {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{}",
            info.filename,
            info.selected_genome,
            info.num_major_variants,
            info.num_minor_variants,
            info.breadth_coverage,
            info.depth_coverage,
            info.num_perfect_kmers,
            info.num_variant_kmers,
            info.num_unmapped_kmers
        ).unwrap();
    }
}


pub fn print_output(
    read_output: &String,
    args: &CallArgs, 
    variants: &Vec<VCFRecord>,
    viral_metadata: &ViralMetadata,
    best_genome_index: &u16, 
){
    info!("Writing output to VCF");

    let file_stem = clean_sample_id(read_output);
        
    let vcf_file = File::create(format!("{}/{}.vcf", args.output, file_stem)).unwrap_or_else(|e| {
        error!("{} | Failed to create vcf output file", e);
        std::process::exit(1);
    });
    let mut writer = BufWriter::new(vcf_file);

    // write out VCF format
    writeln!(writer, "##fileformat=VCFv4.5").unwrap();
    writeln!(writer, "##source=bronko-v{}", BRONKO_VERSION).unwrap();
    writeln!(writer, "##reference=file://{}", read_output).unwrap(); // update to reflect current genome

    let file_meta = &viral_metadata.files[*best_genome_index as usize];

    for seq_entry in &file_meta.sequences {
        let contig = &seq_entry.name;
        let seq_len = seq_entry.len;
        writeln!(writer, "##contig=<ID={},length={}>", contig.split_whitespace().next().unwrap_or(""), seq_len).unwrap();
    }
    writeln!(writer, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">").unwrap();
    writeln!(writer, "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">").unwrap();
    writeln!(writer, "##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Fwd_ref,Rev_ref,Fwd_alt,Rev_alt\">").unwrap();
    writeln!(writer, "##INFO=<ID=SOR,Number=4,Type=Float,Description=\"SOR\">").unwrap();
    writeln!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();

    for variant in variants{
        let seq_out:&str = variant.seq.split_whitespace().next().unwrap_or("");
        writeln!(writer, "{}\t{}\t.\t{}\t{}\t.\tPASS\tDP={};AF={:.3};DP4={},{},{},{};SOR={:.3}", seq_out, variant.pos, nucleotide_bits_to_char(variant.ref_base as u64), nucleotide_bits_to_char(variant.alt_base as u64), variant.depth, variant.af, variant.fwd_ref, variant.rev_ref, variant.fwd_alt, variant.rev_alt, variant.sor).unwrap()
    }
}

#[derive(Debug, Clone)]
pub struct VCFRecord{
    seq: String,
    pos: usize,
    ref_base: u8,
    alt_base: u8,
    fwd_ref: u64,
    rev_ref: u64,
    fwd_alt: u64,
    rev_alt: u64,
    depth: u64,
    af: f64,
    sor: f64
}



pub fn call_variants(
    args: &CallArgs, 
    output: &DashMap<String, OutputData>,
    output_count: &DashMap<String, OutputData>,
    output_rev: &DashMap<String, OutputData>,
    output_rev_count: &DashMap<String, OutputData>,
    min_af: f64,
    filter_end_seq: bool,
    strand_filter: bool,
    n_kmer_per_strand: usize,
    filename: &String, 
) -> (Vec<VCFRecord>, usize, usize, f64, f64) {

    info!("Calling variants for {}", filename);

    let mut results: Vec<VCFRecord> = Vec::new();

    // overall variant calling statitics for sample
    let mut num_minor_variants = 0;
    let mut num_major_variants = 0;

    // overall mapping statistics for sample
    let mut positions_covered: usize = 0; //number of positions with any coverage (used for breadth of coverage)
    let mut total_positions: usize = 0; //length across all sequences
    let mut total_coverage: usize = 0; //total depth across all positions (is averaged out below to get depth of coverage)

    for seq_entry in output.iter(){
        let seq = seq_entry.key();
        let fwd = seq_entry.value();
        let rev = output_rev.get(seq).expect("Missing reverse");
        let fwd_counts = output_count.get(seq).expect("Missing fwd counts");
        let rev_counts = output_rev_count.get(seq).expect("Missing rev counts");

        let len = fwd.counts.len();
        let mut start = 0;
        let mut end = len;
        

        // change the start and end depending if we are filtering the ends of sequences by k
        if filter_end_seq {
            start = args.kmer;
            end = len - args.kmer;
        }

        //add length to total_positions
        total_positions += len as usize;

        //loop through positions
        for i in start..end {
            let row = fwd.counts[i];
            let row_rev = rev.counts[i];
            let count = fwd_counts.counts[i];
            let count_rev = rev_counts.counts[i];


            let ref_base_char = fwd.ref_bases[i];
            let ref_base = nt_to_bits(ref_base_char);

            if ref_base >= 4 {
                continue; // skip non-ACGT
            }

            let pos = i + 1;

            // calculate the depths, including those of fwd and reverse, find the min and max of the two for strand filtering purposes
            let row_total: Vec<u64> = (0..4)
                .map(|b| row[b] + row_rev[b])
                .collect();

            let total_depth = row_total.iter().sum();
            if total_depth == 0 {
                continue; //maybe replace down the line with logic to fix things
            } else {
                positions_covered += 1;
                total_coverage += total_depth as usize;
            }

            // loop through each base and variant call if not the reference base
            for alt_base in 0u8..4 {
                if alt_base == ref_base || row_total[alt_base as usize] == 0 { //skip if the reference base of if there isn't any depth at that position
                    continue;
                }

                // NEW Strand filter logic
                let mut sor = args.strand_odds_max + 1.0;
                if strand_filter {
                    let a = row[ref_base as usize] as f64 + 1.0; //ref fwd
                    let b = row_rev[ref_base as usize] as f64 + 1.0; //ref rev
                    let c = row[alt_base as usize] as f64 + 1.0; //alt fwd
                    let d = row_rev[alt_base as usize] as f64 + 1.0; //alt rev

                    //Using GATK strand odds ratio
                    let r = (a*d)/(b*c);
                    let ref_ratio = (a.min(b)) / (a.max(b));
                    let alt_ratio = (c.min(d)) / (c.max(d));
                    
                    sor = (r + (1.0 / r)).ln() + ref_ratio.ln() - alt_ratio.ln(); 

                    // filter out if greater than strand odds ratio (default 2)
                    if sor > args.strand_odds_max {
                        continue;
                    }

                    // additional filtering for low kmer support across both strands
                    let c_k = count[alt_base as usize] as usize;
                    let d_k = count_rev[alt_base as usize] as usize;

                    if c_k < n_kmer_per_strand && d_k < n_kmer_per_strand {
                        continue;
                    }
                }

                //Get minor af (filter out if below reporting threshold)
                let alt_count = row_total[alt_base as usize];
                let af = alt_count as f64 / total_depth as f64;
                
                if af < min_af {
                    continue;
                }

                // call major/minor variants
                if af >= 0.5 {
                    num_major_variants += 1;
                } else {
                    if total_depth < args.min_depth as u64 { // filter out minor variants when the total depth is too low
                        continue;
                    } 
                    num_minor_variants += 1;
                }

                results.push(VCFRecord { 
                    seq: seq.clone(), 
                    pos: pos, 
                    ref_base: ref_base, 
                    alt_base: alt_base, 
                    fwd_ref: row[ref_base as usize], 
                    rev_ref: row_rev[ref_base as usize],
                    fwd_alt: row[alt_base as usize],
                    rev_alt: row_rev[alt_base as usize],
                    depth: total_depth,
                    af: af,
                    sor: sor
                })

            }
        }

    }

    let breadth_cov: f64 = positions_covered as f64 / total_positions as f64;
    let depth_cov: f64 = total_coverage as f64 / positions_covered as f64;
    info!("Sample breadth of coverage: {}, depth of coverage: {}", breadth_cov, depth_cov);
    info!("Called {} major variants, {} minor above maf = {}", num_major_variants, num_minor_variants, min_af);
    (results, num_major_variants, num_minor_variants, breadth_cov, depth_cov)

}

pub fn count_kmers_kmc(reads: &String, threads: &usize, args: &CallArgs) -> Result<(usize, usize, usize, usize), String> {
    let fastq_path = reads.clone();
    let file_stem = clean_sample_id(&fastq_path);

    let output_dir = Path::new(&args.output);
    let tmp_dir = output_dir.join(format!("tmp_{}", file_stem));
    fs::create_dir_all(&tmp_dir)
        .map_err(|e| format!("Failed to create tmp dir: {}", e))?;

    // make sure parent dir for res_prefix exists
    fs::create_dir_all(output_dir)
        .map_err(|e| format!("Failed to create output dir: {}", e))?;

    let res_prefix = format!("{}/{}.res", args.output, file_stem);
    let kmc_output = Command::new("kmc")
        .args(&[
            &format!("-k{}", args.kmer),
            "-m2",
            &format!("-t{}", threads),
            "-b",
            &format!("-ci{}", args.min_kmers),
            "-cs1000000",
            &format!("{}", fastq_path),
            &res_prefix,
            tmp_dir.to_str().unwrap(),
        ])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .map_err(|e| format!("KMC3 failed: {}", e))?;

    let stdout = String::from_utf8_lossy(&kmc_output.stdout);

    let mut total_reads = None;
    let mut total_kmers = None;
    let mut unique_kmers = None;
    let mut unique_counted_kmers = None;

    for line in stdout.lines() {
        if line.contains("No. of unique counted k-mers") {
            unique_counted_kmers = extract_number_from_line(line);
        } else if line.contains("No. of unique k-mers") {
            unique_kmers = extract_number_from_line(line);
        } else if line.contains("Total no. of k-mers") {
            total_kmers = extract_number_from_line(line);
        } else if line.contains("Total no. of reads") {
            total_reads = extract_number_from_line(line);
        }
    }

    let dump_txt_path = Path::new(&args.output).join(&format!("{}_counts.txt", file_stem));
    let kmc_dump_output = Command::new("kmc_tools")
        .args(&[
            "transform",
            &res_prefix,
            "dump",
            dump_txt_path.to_str().unwrap()
        ])
        .output()
        .map_err(|e| format!("KMC3 tools failed | {}", e))?;

    if let Err(e) = fs::remove_dir_all(&tmp_dir) {
        warn!("Failed to remove tmp dir {}: {}", tmp_dir.display(), e);
    }

    if !kmc_dump_output.status.success(){
        return Err(format!("KMC3 dump failed | {}", String::from_utf8_lossy(&kmc_dump_output.stderr)))
    }

    if let (Some(tr), Some(tk), Some(uk), Some(uck)) = (total_reads, total_kmers, unique_kmers, unique_counted_kmers) {
        Ok((tr, tk, uk, uck))
    } else {
        Err("Failed to parse all KMC stats from stdout".to_string())
    }
}

fn extract_number_from_line(line: &str) -> Option<usize> {
    line.split(':')
        .nth(1)
        .map(|s| s.trim())
        .and_then(|num_str| num_str.replace(',', "").parse::<usize>().ok())
}

#[derive(Clone)]
pub struct OutputData {
    pub counts: Vec<[u64;4]>,
    pub ref_bases: Vec<u8>,
}

fn load_kmers(path: &str) -> Vec<(String, u64)> {
    let file = File::open(path).expect("Failed to open kmer file");
    let reader = BufReader::new(file);

    reader
        .lines()
        .map(|line| {
            let line = line.expect("Read error");
            let mut parts = line.split_whitespace();
            let kmer = parts.next().expect("Missing kmer").to_string();
            let count: u64 = parts.next().expect("Missing count").parse().expect("Invalid count");
            (kmer, count)
        })
        .collect() // Make sure the function's return type is correct!
}

pub fn map_kmers(
    kmers: &Vec<(String, u64)>,
    index: &FxHashMap<u64, Vec<BucketInfo>>,
    viral_metadata: &ViralMetadata,
    threads: &usize,
    args: &CallArgs,
    output_maps: &FxHashMap<
        u16,
        (
            DashMap<String, OutputData>,
            DashMap<String, OutputData>,
            DashMap<String, OutputData>,
            DashMap<String, OutputData>,
        ),
    >,
) -> FxHashMap<u16, (usize, usize, usize)> {
    let k = args.kmer;

    //to store the number of kmers that are mapped perfectly or with 1-edit distance after all of the chunks
    //key is the index of the file in ViralMetadata, values are number of perfectly mapped and 1-edit distance kmers
    let results: DashMap<u16, (usize, usize, usize)> = DashMap::new();

    let chunk_size = if ((kmers.len() / threads) as usize) < 10000 { (kmers.len() / threads) as usize } else { 10000 };

    kmers.par_chunks(chunk_size).for_each(|chunk| {

        //key is the index of the file in ViralMetadata, values are number of perfectly mapped and 1-edit distance kmers and unique perfectly mapped kmers
        let mut local_counts: FxHashMap<u16, (usize, usize, usize)> = FxHashMap::default(); //storing the number of perfectly mapped kmers (buckets found = len(filtered_buckets)) and variant mapped kmers (buckets found = 1) and unique perfect (perfect and only mapped to 1 genome) in thread

        for (kmer, n) in chunk {

            let (kmer_bin, rc) = canonical_kmer(kmer.as_bytes(), k);
            let buckets = assign_buckets(kmer_bin, k);

            let filtered_buckets = if args.use_full_kmer {
                buckets
            } else {
                let len = buckets.len();
                if args.n_fixed * 2 + 1 >= len {
                    vec![] 
                } else {
                    buckets[args.n_fixed..len - args.n_fixed - 1].to_vec()
                }
            };

            let num_buckets_perfect = filtered_buckets.len();
            let mut per_genome_bucket_hits: FxHashMap<u16, usize> = FxHashMap::default(); //number of hits to each genome (where the key is the index of the file) per kmer

            for &bucket in &filtered_buckets {

                if let Some(bucket_infos) = index.get(&bucket) {

                    for info in bucket_infos {
                        // NEED TO UPDATE TO FILTER OUT DUPLICATE BUCKETS IN GENOMES (Problem is that buckets could be from multiple)

                        //get sequence info from metadata
                        let file_meta = &viral_metadata.files[info.file_id as usize];
                        
                        //update the number of hits for this kmer
                        *per_genome_bucket_hits
                            .entry(info.file_id.clone())
                            .or_insert(0) += 1;

                        //get sequence name as well
                        let seq_meta = &file_meta.sequences[info.seq_id as usize];
                        let seq = &seq_meta.name;

                        if let Some(maps) = output_maps.get(&info.file_id) {
                            let (output, output_rev, output_counts, output_rev_counts) = maps;

                            //get genome position and variant position in kmer
                            let genome_pos = info.location as usize;
                            let nuc_x = info.idx as usize;
                            
                            if info.canonical {
                                let pos = k - nuc_x - 1;
                                let bit_idx = (((kmer_bin >> (2 * (k - pos - 1))) & 0b11) ^ 0b11) as usize;
                                let idx  = genome_pos + nuc_x;

                                if rc {
                                    if let Some(mut rec) = output_counts.get_mut(seq) {
                                        rec.counts[idx][bit_idx] += 1;
                                    }
        
                                    if let Some(mut rec) = output.get_mut(seq) {
                                        if rec.counts[idx][bit_idx] < *n {
                                            rec.counts[idx][bit_idx] = *n;
                                        }
                                    }
                                } else {

                                    if let Some(mut rec) = output_rev_counts.get_mut(seq) {
                                        rec.counts[idx][bit_idx] += 1;
                                    }
        
                                    if let Some(mut rec) = output_rev.get_mut(seq) {
                                        if rec.counts[idx][bit_idx] < *n {
                                            rec.counts[idx][bit_idx] = *n;
                                        }
                                    }
                                }
                            } else {
                                let pos = nuc_x;
                                let bit_idx = ((kmer_bin >> (2* (k - pos - 1))) & 0b11) as usize;
                                let idx = genome_pos + nuc_x;

                                if rc {
                                    if let Some(mut rec) = output_rev_counts.get_mut(seq) {
                                        rec.counts[idx][bit_idx] += 1;
                                    }
        
                                    if let Some(mut rec) = output_rev.get_mut(seq) {
                                        if rec.counts[idx][bit_idx] < *n {
                                            rec.counts[idx][bit_idx] = *n;
                                        }
                                    }
                                } else {
                                    if let Some(mut rec) = output_counts.get_mut(seq) {
                                        rec.counts[idx][bit_idx] += 1;
                                    }
        
                                    if let Some(mut rec) = output.get_mut(seq) {
                                        if rec.counts[idx][bit_idx] < *n {
                                            rec.counts[idx][bit_idx] = *n;
                                        }
                                    }
                                }    
                            }
                        }    
                    }
                }
            }

            // Identify perfect + unique hits
            let mut unique_flag: Option<u16> = None;
            for (file_index, hits) in per_genome_bucket_hits.iter() {
                if *hits == num_buckets_perfect {
                    if unique_flag.is_none() {
                        unique_flag = Some(*file_index);
                    } else {
                        // More than one perfect genome → not unique
                        unique_flag = None;
                        break;
                    }
                }
            }

            // Update counts
            for (file_index, hits) in per_genome_bucket_hits {
                let entry = local_counts.entry(file_index).or_insert((0, 0, 0));
                if hits == num_buckets_perfect {
                    entry.0 += 1; // perfect match
                } else if hits > 0 {
                    entry.1 += 1; // variant match
                }
            }

            // If this kmer mapped perfectly to exactly one genome, then update that entry to have one more perfect and unique kmer
            if let Some(unique_genome) = unique_flag {
                let entry = local_counts.entry(unique_genome).or_insert((0, 0, 0));
                entry.2 += 1; // perfect + unique
            }
        }
        // merge local counts into global DashMap
        for (file_index, (perfect, variant, unique_perfect)) in local_counts {
            results
                .entry(file_index)
                .and_modify(|e| {
                    e.0 += perfect;
                    e.1 += variant;
                    e.2 += unique_perfect;
                })
                .or_insert((perfect, variant, unique_perfect));
        }
    });

    results.into_iter().collect()
}


pub fn initialize_output_maps(
    metadata: &ViralMetadata,
) -> FxHashMap<
    u16,
    (
        DashMap<String, OutputData>,
        DashMap<String, OutputData>,
        DashMap<String, OutputData>,
        DashMap<String, OutputData>,
    ),
> {
    let mut result = FxHashMap::default();

    for (i, file) in metadata.files.iter().enumerate() {
        let output = DashMap::new(); //forward depth estimate
        let output_rev = DashMap::new(); //reverse depth estimate
        let output_counts = DashMap::new(); //forward number of kmers
        let output_rev_counts = DashMap::new(); //reverse number of kmers

        for seq_meta in &file.sequences {
            let length = seq_meta.len;
            let mut ref_bases = Vec::with_capacity(length);

            ref_bases.extend_from_slice(&seq_meta.seq);

            let record = OutputData {
                counts: vec![[0u64; 4]; length],
                ref_bases,
            };

            output.insert(seq_meta.name.clone(), record.clone());
            output_rev.insert(seq_meta.name.clone(), record.clone());
            output_counts.insert(seq_meta.name.clone(), record.clone());
            output_rev_counts.insert(seq_meta.name.clone(), record.clone());
        }

        result.insert(
            i as u16,
            (output, output_rev, output_counts, output_rev_counts),
        );
    }

    result
}
