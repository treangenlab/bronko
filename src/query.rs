use crate::cli::*;
use crate::lcb::*;
use crate::util::*;

use anyhow::{Result};
use anyhow::Error;

use needletail::{parse_fastx_file};

use num_cpus;
use log::*;
use crate::consts::*;
use rustc_hash::{FxHashMap, FxHashSet};
use dashmap::DashMap;

use std::fs::File;
use std::sync::{Arc, Mutex};
use std::vec;
use std::fs;
use std::sync::atomic::{AtomicUsize, Ordering};

use rayon::prelude::*;

use std::path::Path;
use std::io::{BufReader, BufRead, BufWriter, Write};
use std::process::{Command, Stdio};

#[derive(Debug, Clone)]
pub struct BucketInfo {
    pub file_name: String, //file name it refers to
    pub seq_name: String, //sequence name (in case mulitple in file)
    pub location: usize, //location of kmer in genome
    pub idx: usize, //nucleotide within the kmer the bucket points to
    pub ref_base: u8, //reference base 00, 01, 10, 11 
    pub canonical: bool, //was the base converted to canonical
}

fn check_args(args: &QueryArgs) {
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

    trace!("{}", args.kmer);

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

    for fasta_file in &args.genomes {
        if !check_fasta(&fasta_file){
            error!("{} does not appear to be a fasta file (must be .fa(.gz)/.fasta(.gz)/.fna(.gz))", &fasta_file);
            std::process::exit(1)
        }
    }

    let available_threads = num_cpus::get();
    if args.threads <= 0 {
        error!("Number of threads must be greater than 0");
        std::process::exit(1);
    } else if args.threads >= available_threads {
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

    // for (r1, r2) in args.first_pairs.iter().zip(args.second_pairs.iter()){
    //     println!("{}, {}", r1, r2)
    // }

}

pub fn query(args: QueryArgs) {

    //First check to make sure none of the arguments are invalid
    check_args(&args);
    trace!("k={}, threads={}", args.kmer, args.threads);

    // create output directory
    let out_path = Path::new(&args.output);
    fs::create_dir(out_path);
    fs::create_dir(out_path.join("tmp"));

    //Get all of the fasta files and process them, get an index
    let (ref_index, seq_info): (FxHashMap<u64, Vec<BucketInfo>>, DashMap<String, usize>) = 
        build_indexes(&args).unwrap_or_else(|e| {
            error!("{} | Reference failed to build", e);
            std::process::exit(1)
        });
    log_memory_usage(true, "Fasta files indexed successfully. Starting counting kmers ");

    // storing the variant information
    let total_samples = args.reads.len() + args.first_pairs.len();
    let mut variant_info: Vec<(String, Vec<VCFRecord>)> = Vec::with_capacity(total_samples); // (sample_name, variant records)

    // PROCESS SINGLE END READS
    if args.reads.len() > 0 {
        for se_read in args.reads.iter() {
            info!("Processing {}", se_read);

            //Get kmer counts
            let (kmers, total_reads, total_kmers, unique_kmers, unique_counted_kmer ) = get_kmers(&se_read, &args);
            info!("{} reads counted from {}", total_reads, se_read);
            info!("{} unique kmers above {} count, {} total unique kmers, {} total kmers (~{} basepairs)", unique_counted_kmer, args.min_kmers, unique_kmers, total_kmers, total_kmers*args.kmer);
            log_memory_usage(true, "Finished counting kmers");

            //initialize output storage and then map the kmers using the index
            let (output, output_count, output_rev, output_rev_count) = initialize_output_maps(&seq_info);
            let (n_variant_mapped, n_perfect_mapped) = map_kmers(&kmers, &ref_index, &args, &output, &output_count, &output_rev, &output_rev_count);
            let message = format!("Mapped {}/{} kmers perfectly, {}/{} had a variant, {} unmapped", n_perfect_mapped, unique_counted_kmer, n_variant_mapped, unique_counted_kmer, unique_counted_kmer-n_perfect_mapped-n_variant_mapped); 
            log_memory_usage(true, &message);
            if ((n_variant_mapped + n_perfect_mapped) as f64 / unique_counted_kmer as f64) < 0.2 {
                warn!("Percent of kmers found is very low, suggesting a bad reference, a bad sequencing run, contamination in sample, or some other issue")
            }

            // call cariants and print them out to vcf
            let variants = call_variants(&args, &output, &output_count, &output_rev, &output_rev_count, args.min_af, !&args.no_end_filter, !args.no_strand_filter, args.n_per_strand);
            log_memory_usage(true, "Called variants successfully");
        
            //print outputs
            if args.output_pileup {
                print_pileup(&se_read, &args, &output, &output_rev, &seq_info);
            }
            print_output(&se_read, &args, &variants, &seq_info);

            variant_info.push((se_read.to_string(), variants));

        }
    }

    // PROCESS PAIRED END READS
    if args.first_pairs.len() > 0 && args.second_pairs.len() > 0 {
        for (r1, r2) in args.first_pairs.iter().zip(args.second_pairs.iter()){
            info!("Processing paired reads {}, {}", r1, r2);
            let (kmers1, total_reads_r1, total_kmers_r1, unique_kmers_r1, unique_counted_kmer_r1 ) = get_kmers(&r1, &args);
            let (kmers2, total_reads_r2, total_kmers_r2, unique_kmers_r2, unique_counted_kmer_r2 ) = get_kmers(&r2, &args);
            info!("{} reads counted from {}", total_reads_r1 + total_reads_r2, r1);
            info!("{} unique kmers above {} count, {} total unique kmers, {} total kmers (~{} basepairs)", unique_counted_kmer_r1 + unique_counted_kmer_r2, args.min_kmers, unique_kmers_r1 + unique_kmers_r2, total_kmers_r1 + total_kmers_r2, (total_kmers_r1 + total_kmers_r2) * args.kmer);
            log_memory_usage(true, "Finished counting kmers");

            //initialize output storage and then map the kmers using the index
            let (output, output_count, output_rev, output_rev_count) = initialize_output_maps(&seq_info);
            let (n_variant_mapped_r1, n_perfect_mapped_r1) = map_kmers(&kmers1, &ref_index, &args, &output, &output_count, &output_rev, &output_rev_count);
            let (n_variant_mapped_r2, n_perfect_mapped_r2) = map_kmers(&kmers2, &ref_index, &args, &output, &output_count, &output_rev, &output_rev_count);
            let message = format!("Mapped {}/{} kmers perfectly, {}/{} had a variant, {} unmapped", n_perfect_mapped_r1 + n_perfect_mapped_r2, unique_counted_kmer_r1 + unique_counted_kmer_r2, n_variant_mapped_r1+n_variant_mapped_r2, unique_counted_kmer_r1 + unique_counted_kmer_r2, (unique_counted_kmer_r1 + unique_counted_kmer_r2)-(n_perfect_mapped_r1+n_perfect_mapped_r2)-(n_variant_mapped_r1 + n_variant_mapped_r2), ); 
            log_memory_usage(true, &message);
            if ((n_variant_mapped_r1 + n_variant_mapped_r2 + n_perfect_mapped_r1 + n_perfect_mapped_r2) as f64 / ((unique_counted_kmer_r1 as f64) + (unique_counted_kmer_r2 as f64))) < 0.2 {
                warn!("Percent of kmers found is very low, suggesting a bad reference, a bad sequencing run, contamination in sample, or some other issue")
            }

            // call cariants and print them out to vcf
            let variants = call_variants(&args, &output, &output_count, &output_rev, &output_rev_count, args.min_af, !&args.no_end_filter, !args.no_strand_filter, args.n_per_strand);
            log_memory_usage(true, "Called variants successfully");

            //print outputs
            if args.output_pileup {
                print_pileup(&r1, &args, &output, &output_rev, &seq_info);
            }
            print_output(&r1, &args, &variants, &seq_info);

            variant_info.push((r1.to_string(), variants));

        }
    }

    info!("All samples processed successfully");

    if args.output_alignment {
        info!("Building alignment");
        build_alignment_fasta(variant_info, &args);
    }

}

pub fn build_alignment_fasta(
    sample_variants: Vec<(String, Vec<VCFRecord>)>,
    args: &QueryArgs
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
    let fasta_out = File::create(format!("{}/alignment.mfa", args.output)).unwrap_or_else(|e| {
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
    writeln!(writer, ">{}", args.genomes[0]).unwrap();
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

pub fn get_kmers(reads_file: &String, args:&QueryArgs) -> (Vec<(String, u64)>, usize, usize, usize, usize){
    //count kmers using kmc3
    let kmc_result = count_kmers_kmc(&reads_file, &args);
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
    args: &QueryArgs,
    output: &DashMap<String, OutputData>,
    output_rev: &DashMap<String, OutputData>,
    seq_info: &DashMap<String, usize>
){
    info!("Writing output to pileup");

    let file_stem = clean_sample_id(read_output);

    let tsv_pileup = File::create(format!("{}/{}.tsv", args.output, file_stem)).unwrap_or_else(|e| {
        error!("{} | Failed to create tsv pileup file", e);
        std::process::exit(1);
    });
    let mut writer = BufWriter::new(tsv_pileup);

    writeln!(writer, "reference\tindex\tref\tA\tC\tG\tT\ta\tc\tg\tt").unwrap();

    for seq_entry in seq_info.iter() {
        let (seq, seq_len) = seq_entry.pair();
        let fwd = output.get(seq).expect("Could not match seq to fwd counts");
        let rev = output_rev.get(seq).expect("Could not match seq to rev counts");

        for (i, (ref_fwd, ref_rev)) in fwd.ref_bases.iter().zip(rev.ref_bases.iter()).enumerate(){
            let ref_base = if *ref_fwd > 3 {*ref_fwd} else {*ref_rev};
            writeln!(writer, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}", seq, &i + 1 , nucleotide_bits_to_char(ref_base as u64), fwd.counts[i][0], fwd.counts[i][1], fwd.counts[i][2], fwd.counts[i][3], rev.counts[i][0], rev.counts[i][1], rev.counts[i][2], rev.counts[i][3]).unwrap();
        }
    }
}


pub fn print_output(
    read_output: &String,
    args: &QueryArgs, 
    variants: &Vec<VCFRecord>,
    seq_info: &DashMap<String, usize>
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
    writeln!(writer, "##source=bronkoV{}", BRONKO_VERSION).unwrap();
    writeln!(writer, "##reference=file://{}", read_output).unwrap(); // update to reflect current genome
    for item in seq_info.iter() {
        let (contig, len) = item.pair();
        writeln!(writer, "##contig=<ID={},length={}>", contig.split_whitespace().next().unwrap_or(""), len).unwrap();
    }
    writeln!(writer, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">").unwrap();
    writeln!(writer, "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">").unwrap();
    writeln!(writer, "##INFO=<ID=DP4,Number=4,Type=Integer,Description=\"Fwd_ref,Rev_ref,Fwd_alt,Rev_alt\">").unwrap();
    writeln!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();

    for variant in variants{
        let seq_out:&str = variant.seq.split_whitespace().next().unwrap_or("");
        writeln!(writer, "{}\t{}\t.\t{}\t{}\t.\tPASS\tDP={};AF={:.3};DP4={},{},{},{}", seq_out, variant.pos, nucleotide_bits_to_char(variant.ref_base as u64), nucleotide_bits_to_char(variant.alt_base as u64), variant.depth, variant.af, variant.fwd_ref, variant.rev_ref, variant.fwd_alt, variant.rev_alt).unwrap()
    }

}

#[derive(Debug)]
struct VCFRecord{
    seq: String,
    pos: usize,
    ref_base: u8,
    alt_base: u8,
    fwd_ref: u64,
    rev_ref: u64,
    fwd_alt: u64,
    rev_alt: u64,
    depth: u64,
    af: f64
}



pub fn call_variants(
    args: &QueryArgs, 
    output: &DashMap<String, OutputData>,
    output_count: &DashMap<String, OutputData>,
    output_rev: &DashMap<String, OutputData>,
    output_rev_count: &DashMap<String, OutputData>,
    min_af: f64,
    filter_end_seq: bool,
    strand_filter: bool,
    n_kmer_per_strand: usize,
) -> Vec<VCFRecord> {

    info!("Calling variants for [replace w filename]");

    let mut results: Vec<VCFRecord> = Vec::new();

    let mut num_minor_variants = 0;
    let mut num_major_variants = 0;

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

        for i in start..end {
            let row = fwd.counts[i];
            let row_rev = rev.counts[i];
            let count = fwd_counts.counts[i];
            let count_rev = rev_counts.counts[i];

            let ref_fwd = fwd.ref_bases[i];
            let ref_rev = rev.ref_bases[i];

            // try to find the ref-base, don't variant call if there is not ACGT in a position
            let ref_base = if ref_fwd < 4 { ref_fwd } else { ref_rev };
            if ref_base >= 4 {
                continue;
            }
            let pos = i + 1;

            // calculate the depths, including those of fwd and reverse, find the min and max of the two for strand filtering purposes
            let row_total: Vec<u64> = (0..4)
                .map(|b| row[b] + row_rev[b])
                .collect();
            let fwd_depth: u64 = row.iter().sum();
            let rev_depth: u64 = row_rev.iter().sum();
            let percent_strand_depth: f64 = args.min_strand_diff;

            let (min_depth_strand, max_depth_strand) = if fwd_depth < rev_depth {
                (fwd_depth, rev_depth)
            } else {
                (rev_depth, fwd_depth)
            };

            let total_depth = row_total.iter().sum();
            if total_depth == 0 {
                continue; //maybe replace down the line with logic to fix things
            }

            // loop through each base and variant call if not the reference base
            for alt_base in 0u8..4 {
                if alt_base == ref_base || row_total[alt_base as usize] == 0 { //skip if the reference base of if there isn't any depth at that position
                    continue;
                }

                /// Strand filter logic
                /// 
                /// If the depths are uneven (one is <min_depth_percent% of the total_depth by default)
                /// then only one of the two strands must pass the n_kmer_per_strand (likely the dominant one)
                /// otherwise both must pass that filter. 
                /// 
                /// If there is no stand filter, then it does not matter, you just let everything pass with the same logic  
                /// 
                let pass_strand_filter = if strand_filter {
                    if min_depth_strand as f64 >= percent_strand_depth * max_depth_strand as f64 {
                        count[alt_base as usize] as usize >= n_kmer_per_strand && count_rev[alt_base as usize] as usize >= n_kmer_per_strand
                    } else {
                        count[alt_base as usize] as usize >= n_kmer_per_strand || count_rev[alt_base as usize] as usize >= n_kmer_per_strand
                    }
                } else {
                    true //might need to change this to follow the portion of above (aka any individual must have n_kmers, but both don't have to)
                };

                if !pass_strand_filter {
                    continue;
                }

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
                    af: af
                })

            }
        }

    }

    info!("Called {} minor + {} major variants above maf={}", num_minor_variants, num_major_variants, min_af);
    results

}

pub fn build_indexes(args: &QueryArgs) -> Result<(FxHashMap<u64, Vec<BucketInfo>>, DashMap<String, usize>), Error> {

    info!("Building indexes from fasta files");
    let k = args.kmer;
    let index: Arc<Mutex<FxHashMap<u64, Vec<BucketInfo>>>> = Arc::new(Mutex::new(FxHashMap::default()));
    let seq_info: DashMap<String, usize> = DashMap::new();

    args.genomes.par_iter().for_each(|file_path|{
        let mut reader = parse_fastx_file(file_path).unwrap_or_else(|e|{
            error!("{} | Failed to parse fasta file: {}", e, file_path);
            std::process::exit(1)
        });

        let file_name = Path::new(file_path).file_stem().and_then(|s| s.to_str()).unwrap_or("unknown").to_string();

        while let Some(record) = reader.next() {
            let record = record.unwrap_or_else(|e|{
                error!("{} | Unable to read record in {}", e, file_path);
                std::process::exit(1)
            });

            let seq: std::borrow::Cow<'_, [u8]> = record.seq();
            let seq_name = String::from_utf8_lossy(record.id()).to_string().split_whitespace().next().unwrap_or_default().to_string();
            let seq_len = seq.len();
            
            seq_info.entry(seq_name.clone()).or_insert(seq_len);

            // trace!("Identified {}", seq_name);
            trace!("{} kmers in {}", seq.len().saturating_sub(k), seq_name);
            for i in 0..=seq.len().saturating_sub(k){
                let kmer: &[u8] = &seq[i..i+k];
                let (kmer_bin, canonical) = canonical_kmer(kmer, k);
                let buckets = assign_buckets(kmer_bin, k);

                for (j, bucket_id) in buckets.iter().enumerate(){
                    let bucket = BucketInfo {
                        file_name: file_name.clone(), //might want to change this in the future as it is kind of redundant
                        seq_name: seq_name.clone(),
                        location: i, //should be able to combine location + idx into single location
                        idx: j,
                        ref_base: nt_to_bits(kmer[j]) as u8, //this also is redundant
                        canonical: canonical //this is maybe the only important thing aside from location
                    };

                    let mut locked = index.lock().unwrap();
                    locked.entry(*bucket_id).or_default().push(bucket);
                }
            }
        }

    });

    Ok((
        Arc::try_unwrap(index)
            .unwrap()
            .into_inner()
            .unwrap(),
        seq_info
    ))
}

pub fn count_kmers_kmc(reads: &String, args: &QueryArgs) -> Result<(usize, usize, usize, usize), String> {
    let fastq_path = reads.clone();
    let file_stem = clean_sample_id(&fastq_path);


    let res_prefix: String= format!("{}/{}.res", args.output, file_stem);
    let kmc_output = Command::new("kmc")
        .args(&[
            &format!("-k{}", args.kmer),
            "-m2",
            &format!("-t{}", args.threads),
            "-b",
            &format!("-ci{}", args.min_kmers),
            "-cs1000000",
            &format!("{}", fastq_path),
            &res_prefix,
            &format!("{}", Path::new(&args.output).join("tmp").to_string_lossy()),
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

    if !kmc_dump_output.status.success(){
        return Err("KMC3 dump failed".to_string())
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
    args: &QueryArgs,
    output: &DashMap<String, OutputData>,
    output_count: &DashMap<String, OutputData>,
    output_rev: &DashMap<String, OutputData>,
    output_rev_count: &DashMap<String, OutputData>
) -> (usize, usize) {
    let k = args.kmer;

    let num_variant_mapped = Arc::new(AtomicUsize::new(0));
    let num_perfect_mapped = Arc::new(AtomicUsize::new(0));

    let variant_mapped = Arc::clone(&num_variant_mapped);
    let perfect_mapped = Arc::clone(&num_perfect_mapped);

    let chunk_size = if ((kmers.len() / args.threads) as usize) < 10000 { (kmers.len() / args.threads) as usize } else { 10000 };

    kmers.par_chunks(chunk_size).for_each(|chunk| {

        for (kmer, n) in chunk {

            let (kmer_bin, rc) = canonical_kmer(kmer.as_bytes(), k);
            let buckets = assign_buckets(kmer_bin, k);

            let mut found_buckets = 0;

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

            for &bucket in &filtered_buckets {

                if let Some(bucket_infos) = index.get(&bucket) {
                    found_buckets += 1;

                    if bucket_infos.len() == 1 {
                        let info = &bucket_infos[0];
                        let genome_pos = info.location;
                        let seq = &info.seq_name;
                        let nuc_x = info.idx;
                        let ref_base = info.ref_base;
                        
                        if info.canonical {
                            let pos = k - nuc_x - 1;
                            let bit_idx = (((kmer_bin >> (2 * (k - pos - 1))) & 0b11) ^ 0b11) as usize;
                            let idx  = genome_pos + nuc_x;

                            if rc {
                                if let Some(mut rec) = output_count.get_mut(seq) {
                                    rec.counts[idx][bit_idx] += 1;
                                }
    
                                if let Some(mut rec) = output.get_mut(seq) {
                                    if rec.counts[idx][bit_idx] < *n {
                                        rec.counts[idx][bit_idx] = *n;
                                        rec.ref_bases[idx] = ref_base;
                                    }
                                }
                            } else {

                                if let Some(mut rec) = output_rev_count.get_mut(seq) {
                                    rec.counts[idx][bit_idx] += 1;
                                }
    
                                if let Some(mut rec) = output_rev.get_mut(seq) {
                                    if rec.counts[idx][bit_idx] < *n {
                                        rec.counts[idx][bit_idx] = *n;
                                        rec.ref_bases[idx] = ref_base;
                                    }
                                }
                            }
                        } else {
                            let pos = nuc_x;
                            let bit_idx = ((kmer_bin >> (2* (k - pos - 1))) & 0b11) as usize;
                            let idx = genome_pos + nuc_x;

                            if rc {
                                if let Some(mut rec) = output_rev_count.get_mut(seq) {
                                    rec.counts[idx][bit_idx] += 1;
                                }
    
                                if let Some(mut rec) = output_rev.get_mut(seq) {
                                    if rec.counts[idx][bit_idx] < *n {
                                        rec.counts[idx][bit_idx] = *n;
                                        rec.ref_bases[idx] = ref_base;
                                    }
                                }
                            } else {
                                if let Some(mut rec) = output_count.get_mut(seq) {
                                    rec.counts[idx][bit_idx] += 1;
                                }
    
                                if let Some(mut rec) = output.get_mut(seq) {
                                    if rec.counts[idx][bit_idx] < *n {
                                        rec.counts[idx][bit_idx] = *n;
                                        rec.ref_bases[idx] = ref_base;
                                    }
                                }
                            }
                        }
                        

                    }
                }
            }

            if found_buckets == 1 {
                variant_mapped.fetch_add(1, Ordering::Relaxed);
            }
    
            if found_buckets == k {
                perfect_mapped.fetch_add(1, Ordering::Relaxed);
            }
        }
    });

    (
        num_variant_mapped.load(Ordering::Relaxed),
        num_perfect_mapped.load(Ordering::Relaxed),
    )
}


pub fn initialize_output_maps(
    seq_info: &DashMap<String, usize>,
) -> (
    DashMap<String, OutputData>,
    DashMap<String, OutputData>,
    DashMap<String, OutputData>,
    DashMap<String, OutputData>,
) {
    let output = DashMap::new();
    let output_rev = DashMap::new();
    let output_counts = DashMap::new();
    let output_rev_counts = DashMap::new();

    for entry in seq_info.iter() {
        let seq_name = entry.key().clone();
        let length = *entry.value();

        let record = OutputData {
            counts: vec![[0u64; 4]; length],  // vector of 4 counts per base position
            ref_bases: vec![b' '; length],   // prefill with spaces, will hold the reference at that position
        };

        output.insert(seq_name.clone(), record.clone());
        output_rev.insert(seq_name.clone(), record.clone());
        output_counts.insert(seq_name.clone(), record.clone());
        output_rev_counts.insert(seq_name.clone(), record.clone());
    }

    (output, output_rev, output_counts, output_rev_counts)
}