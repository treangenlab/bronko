use crate::cli::*;
use crate::lcb::*;
use crate::util::*;
use crate::build::*;
use crate::indels::*;
use crate::consts::*;

use anyhow::{Result};

use num_cpus;
use log::*;
use rustc_hash::{FxHashMap, FxHashSet};
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
use std::cmp::max;

use statrs::distribution::{StudentsT, ContinuousCDF};



fn check_args(args: &CallArgs) {
    let output_level;
    if args.verbose {
        output_level = log::LevelFilter::Trace;
    } else if args.debug {
        output_level = log::LevelFilter::Debug;
    } else {
        output_level = log::LevelFilter::Info;
    }

    // Log file lives inside the run's output folder and mirrors the terminal.
    let log_path = Path::new(&args.output).join("bronko.log");
    init_logging(output_level, &log_path);
    print_banner();

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

    if let Some(file_input_path) = &args.file_input {
        if check_txt(file_input_path){
            if !Path::new(file_input_path).exists(){
                error!("Filepath provided with --file-input does not exist: {}", file_input_path);
                std::process::exit(1);
            }
        } else {
            error!("Filepath provided with --file-input does not appear to be a .txt file: {}", file_input_path);
            std::process::exit(1);
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


    if args.strand_balance_ratio < 0.0 {
        error!("Strand balance ratio is set to below 0, must be between 0.0 and 1.0");
        std::process::exit(1);
    } else if args.strand_balance_ratio > 1.0 {
        error!("Strand balance ratio is set above 1, must be between 0.0 and 1.0");
        std::process::exit(1);
    } else if args.strand_balance_ratio == 1.0 {
        warn!("Strand balance ratio is set to 1, all variants will pass this filter");
    }

    if args.min_variant_depth < 0 {
        warn!("Minimum variant depth set below 0, all variants will be returned if passing other thresholds");
    }

    if args.min_depth < 0 {
        warn!("Minimum total depth for minor variant calling set below 0, all variants will be returned if passing other thresholds");
    }

    if args.variant_multiplier < 1.0 {
        error!("Noise multiplier for variant detection is set to below 1.0, must be greater than 1.0 (recommended between 1.3-2.0)");
        std::process::exit(1);
    } else if args.variant_multiplier > 2.0 {
        warn!("Strand balance ratio is set above 2, may experience a drop in recall (we recommend ~1.5)");
    } else if args.variant_multiplier == 1.0 {
        warn!("Noise multiplier for variant detection set to 1.0, all variants will pass this filter");
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
    num_major_snp: usize,
    num_minor_snp: usize,
    num_insertions: usize,
    num_deletions: usize,
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

    //collect all the fastq paths from the file and other sources 
    let mut n_se_seq_datasets = args.reads.len();
    let mut n_pe_seq_datasets = args.first_pairs.len();
    let mut se_seq_datasets = if n_se_seq_datasets == 0 {
        Vec::new()
    } else {
        canonicalize_file_paths(&args.reads)
    };
    
    let mut pe_seq_datasets1 = if n_pe_seq_datasets == 0 {
        Vec::new()
    } else {
        canonicalize_file_paths(&args.first_pairs)
    };

    let mut pe_seq_datasets2 = if n_pe_seq_datasets == 0 {
        Vec::new()
    } else {
        canonicalize_file_paths(&args.second_pairs)
    };
    
    if n_pe_seq_datasets == 0 {
        trace!("No PE datasets provided with -1/-2");
    }

    if n_se_seq_datasets ==  0 {
        trace!("No SE datasets provided with -r");
    }

    if let Some(file_input_path) = &args.file_input {
        let genomes_file = File::open(file_input_path).unwrap_or_else(|e| {
            error!("{} | Failed to open genomes file", e);
            std::process::exit(1);
        });

        let reader = BufReader::new(genomes_file);
                
        for line in reader.lines() {
            match line {
                Ok(line_path) => {
                    //split by tab to get if paired or not
                    let parts: Vec<&str> = line_path.split_whitespace().collect();

                    match parts.as_slice() {
                        //SE reads
                        [path] => {
                            let p = path.to_string();
                            check_fastq(&p);

                            let canonical_path = fs::canonicalize(&p).unwrap_or_else(|e| {
                                error!("{} | Path within txt file does not exist or is inaccessible: {}", e, p);
                                std::process::exit(1);
                            });
                            let path_string = canonical_path.to_string_lossy().into_owned();

                            if !se_seq_datasets.contains(&path_string) {
                                se_seq_datasets.push(path_string);
                            }
                            n_se_seq_datasets += 1;
                        }

                        //PE reads
                        [path1, path2] => {
                            let p1 = path1.to_string();
                            let p2 = path2.to_string();

                            check_fastq(&p1);
                            check_fastq(&p2);

                            let canonical_path1 = fs::canonicalize(&p1).unwrap_or_else(|e| {
                                error!("{} | Path within txt file does not exist or is inaccessible: {}", e, p1);
                                std::process::exit(1);
                            });
                            let path_string1 = canonical_path1.to_string_lossy().into_owned();

                            let canonical_path2 = fs::canonicalize(&p2).unwrap_or_else(|e| {
                                error!("{} | Path within txt file does not exist or is inaccessible: {}", e, p2);
                                std::process::exit(1);
                            });
                            let path_string2 = canonical_path2.to_string_lossy().into_owned();

                            if !pe_seq_datasets1.contains(&path_string1) && 
                               !pe_seq_datasets2.contains(&path_string1) && 
                               !pe_seq_datasets1.contains(&path_string2) &&
                               !pe_seq_datasets2.contains(&path_string2) 
                            {
                                pe_seq_datasets1.push(path_string1);
                                pe_seq_datasets2.push(path_string2);
                                n_pe_seq_datasets += 1;
                            }
                        }

                        _ => {
                            error!("Invalid format in genome list: {}", line_path);
                            std::process::exit(1);
                        }
                    }
                }
                Err(e) => {
                    error!("{} | Failed to read line in {}", e, file_input_path);
                    std::process::exit(1);
                }
            }
        }
    }

    if n_pe_seq_datasets + n_se_seq_datasets == 0 {
        error!("No sequencing datasets provided using -r or -1/-2 or --file-input, exiting");
        std::process::exit(1);
    } else {
        info!("Read in {} sequencing datasets ({} single end, {} paired end)", n_pe_seq_datasets + n_se_seq_datasets, n_se_seq_datasets, n_pe_seq_datasets);
    }

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
    let total_samples = n_se_seq_datasets + n_pe_seq_datasets;
    let mut variant_info: Vec<(String, Vec<VCFRecord>)> = Vec::with_capacity(total_samples); // (sample_name, variant records)
    let mut output_info: Vec<OutputInfo> = Vec::with_capacity(total_samples); //storing the outputs for tsv overview

    // PROCESS SINGLE END READS
    if n_se_seq_datasets > 0 {
        for se_read in se_seq_datasets.iter() {
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
            let flanking_base_map: DashMap<u64, [[u64; 4]; 2]> = DashMap::new();
            let mapping_data = map_kmers(&kmers, &ref_index, &viral_metadata, &args.threads, &args, &output_maps, &flanking_base_map);
            trace!("Collected flanking base counts for {} buckets", flanking_base_map.len());

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

            
            // first pass to call variants
            let (output, output_rev, output_counts, output_rev_counts)= output_maps.get(&best_genome_index).unwrap_or_else(|| {
                error!("Failed to find mapping data for selected genome");
                std::process::exit(1);
            });

            let (first_pass_variants, fp_major, fp_minor, fp_breadth, fp_depth) = call_variants(
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
                None, // first pass: original-reference coordinates
            );


            // reconstruct indels and sequence ends unless disabled
            let (reconstructed_indels, reconstructed_ends) = if args.no_indels {
                info!("Indel reconstruction disabled (--no-indels)");
                (Vec::new(), Vec::new())
            } else {
                let infos = identify_indel_breakpoints(output, output_rev, args.indel_max_ratio, args.indel_min_depth as u64, args.indel_max_drop_depth as u64, args.indel_width);
                let (merged_pairs, end_anchors) = plan_indel_events(&infos, args.kmer, args.indel_min_depth as u64);
                let n_breakpoints: usize = infos.iter().map(|i| i.breakpoints.len()).sum();
                let n_pairs: usize = merged_pairs.iter().map(|(_, p)| p.len()).sum();
                info!("Identified {} breakpoints across {} sequences; planned {} internal indel event(s) and {} end handoff(s)", n_breakpoints, infos.len(), n_pairs, end_anchors.len());
                trace!("Planned internal indel pairs: {:?}", merged_pairs);
                trace!("Planned end anchors: {:?}", end_anchors);

                let reconstructed_indels = reconstruct_indels(output, &flanking_base_map, &merged_pairs, args.kmer);
                info!("Reconstructed {} indel alleles", reconstructed_indels.len());
                for indel in &reconstructed_indels {
                    info!("  indel {}:{}-{} (ref {}-{}) allele: {}", indel.seq, indel.drop_pos, indel.rise_pos, indel.ref_start, indel.ref_end, String::from_utf8_lossy(&indel.allele));
                }

                let reconstructed_ends = reconstruct_ends(output, output_rev, &flanking_base_map, args.kmer, args.indel_min_depth as u64, &end_anchors);
                info!("Reconstructed {} sequence end(s) for the consensus", reconstructed_ends.len());
                for end in &reconstructed_ends {
                    info!("  end {} (ref {}-{}) allele: {}", end.seq, end.ref_start, end.ref_end, String::from_utf8_lossy(&end.allele));
                }

                (reconstructed_indels, reconstructed_ends)
            };

            // rebuild a consensus from the indels and ends
            let consensus_edits: Vec<ReconstructedIndel> = reconstructed_indels
                .iter()
                .cloned()
                .chain(reconstructed_ends.iter().filter(|end| {
                    !reconstructed_indels.iter().any(|ind| ind.seq == end.seq
                        && end.ref_start <= ind.ref_end && ind.ref_start <= end.ref_end)
                }).cloned())
                .collect();

            let has_indels = !consensus_edits.is_empty();

            if !has_indels {
                info!("No indels or sequence ends were reconstructed, skipping second pass of mapping");
            } else {
                info!("Printing reconstructed consensus sequence for sample {} ({} indel/ends)", clean_sample_id(se_read), consensus_edits.len());
            }

            // print the sample consensus so we can use it for the second pass (if anything was reconstructed)
            if args.output_consensus || has_indels {
                print_consensus(&se_read, &args, &output, &output_rev, &viral_metadata, &best_genome_index, &consensus_edits, &first_pass_variants);
            }

            // Second pass of mapping if reconstructed
            let second_pass = if has_indels {
                let consensus_path = format!("{}/{}.fa", args.output, clean_sample_id(se_read));
                let (new_index, new_meta) = build_indexes(args.kmer, &[consensus_path]).unwrap_or_else(|e| {
                    error!("{} | Failed to build index from corrected consensus", e);
                    std::process::exit(1);
                });
                info!("Built corrected index from consensus: {} buckets, {} sequence(s)", new_index.len(), new_meta.files.iter().map(|f| f.sequences.len()).sum::<usize>());

                // flanking base counts aren't needed here, so the pass-1 map is reused as a throwaway sink.
                let new_output_maps = initialize_output_maps(&new_meta);
                let new_mapping_data = map_kmers(&kmers, &new_index, &new_meta, &args.threads, &args, &new_output_maps, &flanking_base_map);
                let (new_perfect, new_variant, _) = new_mapping_data.get(&0).copied().unwrap_or((0, 0, 0));
                info!("Re-mapped kmers onto corrected consensus: {}/{} perfect, {} variant", new_perfect, unique_counted_kmer, new_variant);

                // map each corrected coordinate back to the original reference (positions inside a spliced
                // indel allele have no original coordinate and are skipped during calling)
                let mut coord_maps: FxHashMap<String, CoordMap> = FxHashMap::default();
                for entry in output.iter() {
                    let seq = entry.key();
                    let ref_len = entry.value().ref_bases.len();
                    coord_maps.insert(seq.clone(), build_coord_map(ref_len, seq, &reconstructed_indels));
                }

                Some((new_meta, new_output_maps, coord_maps, (new_perfect, new_variant)))
            } else {
                None
            };

            // pick the pileup source + metadata all output runs against: the corrected second pass when
            // present, otherwise the selected genome's first-pass pileup
            let (cur_output, cur_output_rev, cur_output_counts, cur_output_rev_counts, cur_meta, cur_genome_index, coord_maps_opt) =
                match &second_pass {
                    Some((new_meta, new_output_maps, coord_maps, _)) => {
                        let (o, orv, oc, orc) = new_output_maps.get(&0).unwrap_or_else(|| {
                            error!("Failed to find mapping data for corrected consensus");
                            std::process::exit(1);
                        });
                        (o, orv, oc, orc, new_meta, 0u16, Some(coord_maps))
                    }
                    None => (output, output_rev, output_counts, output_rev_counts, &viral_metadata, best_genome_index, None),
                };

            // Print the final variants. First-pass major SNVs are always carried over. If there was a second pass, the
            // the minor variants and any indels or new SNVs are also carried over and translated position-wise. Without a second pass, the first-pass
            // call is the final result as-is.
            let (mut variants, num_major_snp, num_minor_snp, breadth_cov, depth_cov) = if let Some(coord_maps) = coord_maps_opt {
                let (sp_variants, _sp_major, _sp_minor, sp_breadth, sp_depth) = call_variants(
                    &args,
                    cur_output,
                    cur_output_counts,
                    cur_output_rev,
                    cur_output_rev_counts,
                    args.min_af,
                    !args.no_end_filter,
                    !args.no_strand_filter,
                    args.n_per_strand,
                    &best_genome_filename,
                    coord_maps_opt, // second pass: translate corrected coords -> original reference
                );

                let mut merged: Vec<VCFRecord> = Vec::new();
                let mut major_keys: FxHashSet<(String, usize)> = FxHashSet::default();
                for r in &first_pass_variants {
                    if !(r.af >= 0.5 && r.ref_allele.len() == 1 && r.alt_allele.len() == 1) {
                        continue; // only major SNVs are carried over from the first pass
                    }
                    major_keys.insert((r.seq.clone(), r.pos));
                    let cm = match coord_maps.get(&r.seq) {
                        Some(c) => c,
                        None => continue,
                    };
                    // a position with no corrected coordinate lies inside a reconstructed indel and is
                    // subsumed by that indel's record
                    let corrected = match cm.orig_to_corrected.get(r.pos - 1).copied().flatten() {
                        Some(c) => c,
                        None => continue,
                    };
                    // recompute against the second-pass pileup; if it has no coverage there, keep the
                    // first-pass call
                    match major_snv_second_pass_record(&r.seq, r.pos, corrected, r.ref_allele[0], r.alt_allele[0], cur_output, cur_output_rev) {
                        Some(rec) => merged.push(rec),
                        None => merged.push(r.clone()),
                    }
                }
                // add the second-pass calls except where a first-pass major already reports that position
                for r in sp_variants {
                    if !major_keys.contains(&(r.seq.clone(), r.pos)) {
                        merged.push(r);
                    }
                }

                // recount SNVs from the merged set (indel records are added and counted separately below)
                let num_major_snp = merged.iter().filter(|r| r.af >= 0.5 && r.ref_allele.len() == 1 && r.alt_allele.len() == 1).count();
                let num_minor_snp = merged.iter().filter(|r| r.af < 0.5 && r.ref_allele.len() == 1 && r.alt_allele.len() == 1).count();

                (merged, num_major_snp, num_minor_snp, sp_breadth, sp_depth)
            } else {
                (first_pass_variants, fp_major, fp_minor, fp_breadth, fp_depth)
            };
            log_memory_usage(true, "Called variants successfully");

            // emit each reconstructed indel as its own VCF record (original-reference coordinates, depth
            // pulled from the corrected-consensus pileup), tallying insertions and deletions separately
            let mut num_insertions = 0;
            let mut num_deletions = 0;
            if let Some(coord_maps) = coord_maps_opt {
                for indel in &reconstructed_indels {
                    if let Some(coord_map) = coord_maps.get(&indel.seq) {
                        if let Some(rec) = indel_to_vcf_record(indel, output, output_rev, cur_output, cur_output_rev, coord_map) {
                            match rec.var_type() {
                                "INS" => num_insertions += 1,
                                "DEL" => num_deletions += 1,
                                _ => {}
                            }
                            variants.push(rec);
                        }
                    }
                }
            }

            //print outputs
            if args.output_pileup {
                print_pileup(&se_read, &args, cur_output, cur_output_rev, cur_meta, &cur_genome_index, second_pass.is_some());
            }

            // VCF reports in original-reference coordinates regardless of pass
            print_output(&se_read, &args, &variants, &viral_metadata, &best_genome_index);

            // kmer mapping statistics: use the corrected second pass when indels were applied, else the
            // first pass against the selected genome
            let (out_perfect_kmers, out_variant_kmers, out_unmapped_kmers) = match &second_pass {
                Some((.., (p, v))) => (*p, *v, unique_counted_kmer.saturating_sub(*p + *v)),
                None => (*n_perfect_mapped, *n_variant_mapped, n_unmapped_kmers),
            };

            output_info.push(OutputInfo {
                filename: se_read.to_string(),
                selected_genome: best_genome_filename.to_string(),
                num_major_snp: num_major_snp,
                num_minor_snp: num_minor_snp,
                num_insertions: num_insertions,
                num_deletions: num_deletions,
                breadth_coverage: breadth_cov,
                depth_coverage: depth_cov,
                num_perfect_kmers: out_perfect_kmers,
                num_variant_kmers: out_variant_kmers,
                num_unmapped_kmers: out_unmapped_kmers,
            });
            variant_info.push((se_read.to_string(), variants));

            cleanup_kmc_files(&args.output, &args.keep_kmer_counts);

        }
    }

    // PROCESS PAIRED END READS
    if n_pe_seq_datasets > 0 {
        for (r1, r2) in pe_seq_datasets1.iter().zip(pe_seq_datasets2.iter()){
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
            let flanking_base_map: DashMap<u64, [[u64; 4]; 2]> = DashMap::new();
            let mapping_data_r1 = map_kmers(&kmers1, &ref_index, &viral_metadata, &half_threads, &args, &output_maps, &flanking_base_map);
            let mapping_data_r2 = map_kmers(&kmers2, &ref_index, &viral_metadata, &half_threads, &args, &output_maps, &flanking_base_map);
            trace!("Collected flanking base counts for {} buckets", flanking_base_map.len());

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

            // call variants for first pass
            let (output, output_rev, output_counts, output_rev_counts)= output_maps.get(&best_genome_index).unwrap_or_else(|| {
                error!("Failed to find mapping data for selected genome");
                std::process::exit(1);
            });

            let (first_pass_variants, fp_major, fp_minor, fp_breadth, fp_depth) = call_variants(
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
                None, // first pass: original-reference coordinates
            );


            // reconstruct indels and sequence ends unless disabled
            let (reconstructed_indels, reconstructed_ends) = if args.no_indels {
                info!("Indel reconstruction disabled (--no-indels)");
                (Vec::new(), Vec::new())
            } else {
                let infos = identify_indel_breakpoints(output, output_rev, args.indel_max_ratio, args.indel_min_depth as u64, args.indel_max_drop_depth as u64, args.indel_width);
                let (merged_pairs, end_anchors) = plan_indel_events(&infos, args.kmer, args.indel_min_depth as u64);
                let n_breakpoints: usize = infos.iter().map(|i| i.breakpoints.len()).sum();
                let n_pairs: usize = merged_pairs.iter().map(|(_, p)| p.len()).sum();
                info!("Identified {} breakpoints across {} sequences; planned {} internal indel event(s) and {} end handoff(s)", n_breakpoints, infos.len(), n_pairs, end_anchors.len());
                trace!("Planned internal indel pairs: {:?}", merged_pairs);
                trace!("Planned end anchors: {:?}", end_anchors);

                let reconstructed_indels = reconstruct_indels(output, &flanking_base_map, &merged_pairs, args.kmer);
                info!("Reconstructed {} indel alleles", reconstructed_indels.len());
                for indel in &reconstructed_indels {
                    info!("  indel {}:{}-{} (ref {}-{}) allele: {}", indel.seq, indel.drop_pos, indel.rise_pos, indel.ref_start, indel.ref_end, String::from_utf8_lossy(&indel.allele));
                }

                let reconstructed_ends = reconstruct_ends(output, output_rev, &flanking_base_map, args.kmer, args.indel_min_depth as u64, &end_anchors);
                info!("Reconstructed {} sequence end(s) for the consensus", reconstructed_ends.len());
                for end in &reconstructed_ends {
                    info!("  end {} (ref {}-{}) allele: {}", end.seq, end.ref_start, end.ref_end, String::from_utf8_lossy(&end.allele));
                }

                (reconstructed_indels, reconstructed_ends)
            };

            // rebuild a consensus from the indels and ends
            let consensus_edits: Vec<ReconstructedIndel> = reconstructed_indels
                .iter()
                .cloned()
                .chain(reconstructed_ends.iter().filter(|end| {
                    !reconstructed_indels.iter().any(|ind| ind.seq == end.seq
                        && end.ref_start <= ind.ref_end && ind.ref_start <= end.ref_end)
                }).cloned())
                .collect();

            let has_indels = !consensus_edits.is_empty();

            if !has_indels {
                info!("No indels or sequence ends were reconstructed, skipping second pass of mapping");
            } else {
                info!("Printing reconstructed consensus sequence for sample {} ({} indel/ends)", clean_sample_id(r1), consensus_edits.len());
            }

            // print the sample consensus so we can use it for the second pass (if anything was reconstructed)
            if args.output_consensus || has_indels {
                print_consensus(&r1, &args, &output, &output_rev, &viral_metadata, &best_genome_index, &consensus_edits, &first_pass_variants);
            }

            // Second pass of mapping if reconstructed
            let second_pass = if has_indels {
                let consensus_path = format!("{}/{}.fa", args.output, clean_sample_id(r1));
                let (new_index, new_meta) = build_indexes(args.kmer, &[consensus_path]).unwrap_or_else(|e| {
                    error!("{} | Failed to build index from corrected consensus", e);
                    std::process::exit(1);
                });
                info!("Built corrected index from consensus: {} buckets, {} sequence(s)", new_index.len(), new_meta.files.iter().map(|f| f.sequences.len()).sum::<usize>());

                // flanking base counts aren't needed here, so the pass-1 map is reused as a throwaway sink.
                let new_output_maps = initialize_output_maps(&new_meta);
                let new_mapping_data_r1 = map_kmers(&kmers1, &new_index, &new_meta, &half_threads, &args, &new_output_maps, &flanking_base_map);
                let new_mapping_data_r2 = map_kmers(&kmers2, &new_index, &new_meta, &half_threads, &args, &new_output_maps, &flanking_base_map);
                let (p1, v1, _) = new_mapping_data_r1.get(&0).copied().unwrap_or((0, 0, 0));
                let (p2, v2, _) = new_mapping_data_r2.get(&0).copied().unwrap_or((0, 0, 0));
                let (new_perfect, new_variant) = (p1 + p2, v1 + v2);
                info!("Re-mapped kmers onto corrected consensus: {}/{} perfect, {} variant", new_perfect, unique_counted_kmer_r1 + unique_counted_kmer_r2, new_variant);

                // map each corrected coordinate back to the original reference (positions inside a spliced
                // indel allele have no original coordinate and are skipped during calling)
                let mut coord_maps: FxHashMap<String, CoordMap> = FxHashMap::default();
                for entry in output.iter() {
                    let seq = entry.key();
                    let ref_len = entry.value().ref_bases.len();
                    coord_maps.insert(seq.clone(), build_coord_map(ref_len, seq, &reconstructed_indels));
                }

                Some((new_meta, new_output_maps, coord_maps, (new_perfect, new_variant)))
            } else {
                None
            };

            // pick the pileup source + metadata all output runs against: the corrected second pass when
            // present, otherwise the selected genome's first-pass pileup
            let (cur_output, cur_output_rev, cur_output_counts, cur_output_rev_counts, cur_meta, cur_genome_index, coord_maps_opt) =
                match &second_pass {
                    Some((new_meta, new_output_maps, coord_maps, _)) => {
                        let (o, orv, oc, orc) = new_output_maps.get(&0).unwrap_or_else(|| {
                            error!("Failed to find mapping data for corrected consensus");
                            std::process::exit(1);
                        });
                        (o, orv, oc, orc, new_meta, 0u16, Some(coord_maps))
                    }
                    None => (output, output_rev, output_counts, output_rev_counts, &viral_metadata, best_genome_index, None),
                };
            
            // Print the final variants. First-pass major SNVs are always carried over. If there was a second pass, the
            // the minor variants and any indels or new SNVs are also carried over and translated position-wise. Without a second pass, the first-pass
            // call is the final result as-is.
            let (mut variants, num_major_snp, num_minor_snp, breadth_cov, depth_cov) = if let Some(coord_maps) = coord_maps_opt {
                let (sp_variants, _sp_major, _sp_minor, sp_breadth, sp_depth) = call_variants(
                    &args,
                    cur_output,
                    cur_output_counts,
                    cur_output_rev,
                    cur_output_rev_counts,
                    args.min_af,
                    !args.no_end_filter,
                    !args.no_strand_filter,
                    args.n_per_strand,
                    &best_genome_filename,
                    coord_maps_opt, // second pass: translate corrected coords -> original reference
                );

                let mut merged: Vec<VCFRecord> = Vec::new();
                let mut major_keys: FxHashSet<(String, usize)> = FxHashSet::default();
                for r in &first_pass_variants {
                    if !(r.af >= 0.5 && r.ref_allele.len() == 1 && r.alt_allele.len() == 1) {
                        continue; // only major SNVs are carried over from the first pass
                    }
                    major_keys.insert((r.seq.clone(), r.pos));
                    let cm = match coord_maps.get(&r.seq) {
                        Some(c) => c,
                        None => continue,
                    };
                    // a position with no corrected coordinate lies inside a reconstructed indel and is
                    // subsumed by that indel's record
                    let corrected = match cm.orig_to_corrected.get(r.pos - 1).copied().flatten() {
                        Some(c) => c,
                        None => continue,
                    };
                    // recompute against the second-pass pileup; if it has no coverage there, keep the
                    // first-pass call
                    match major_snv_second_pass_record(&r.seq, r.pos, corrected, r.ref_allele[0], r.alt_allele[0], cur_output, cur_output_rev) {
                        Some(rec) => merged.push(rec),
                        None => merged.push(r.clone()),
                    }
                }
                // add the second-pass calls except where a first-pass major already reports that position
                for r in sp_variants {
                    if !major_keys.contains(&(r.seq.clone(), r.pos)) {
                        merged.push(r);
                    }
                }

                // recount SNVs from the merged set (indel records are added and counted separately below)
                let num_major_snp = merged.iter().filter(|r| r.af >= 0.5 && r.ref_allele.len() == 1 && r.alt_allele.len() == 1).count();
                let num_minor_snp = merged.iter().filter(|r| r.af < 0.5 && r.ref_allele.len() == 1 && r.alt_allele.len() == 1).count();

                (merged, num_major_snp, num_minor_snp, sp_breadth, sp_depth)
            } else {
                (first_pass_variants, fp_major, fp_minor, fp_breadth, fp_depth)
            };
            log_memory_usage(true, "Called variants successfully");

            // emit each reconstructed indel as its own VCF record (original-reference coordinates, depth
            // pulled from the corrected-consensus pileup), tallying insertions and deletions separately
            let mut num_insertions = 0;
            let mut num_deletions = 0;
            if let Some(coord_maps) = coord_maps_opt {
                for indel in &reconstructed_indels {
                    if let Some(coord_map) = coord_maps.get(&indel.seq) {
                        if let Some(rec) = indel_to_vcf_record(indel, output, output_rev, cur_output, cur_output_rev, coord_map) {
                            match rec.var_type() {
                                "INS" => num_insertions += 1,
                                "DEL" => num_deletions += 1,
                                _ => {}
                            }
                            variants.push(rec);
                        }
                    }
                }
            }

            //print outputs
            if args.output_pileup {
                print_pileup(&r1, &args, cur_output, cur_output_rev, cur_meta, &cur_genome_index, second_pass.is_some());
            }

            // VCF reports in original-reference coordinates regardless of pass
            print_output(&r1, &args, &variants, &viral_metadata, &best_genome_index);

            // kmer mapping statistics: use the corrected second pass when indels were applied, else the
            // first pass against the selected genome
            let (out_perfect_kmers, out_variant_kmers, out_unmapped_kmers) = match &second_pass {
                Some((.., (p, v))) => (*p, *v, (unique_counted_kmer_r1 + unique_counted_kmer_r2).saturating_sub(*p + *v)),
                None => (*n_perfect_mapped_r1 + *n_perfect_mapped_r2, *n_variant_mapped_r1 + *n_variant_mapped_r2, n_unmapped_kmers),
            };

            output_info.push(OutputInfo {
                filename: r1.to_string(),
                selected_genome: best_genome_filename.to_string(),
                num_major_snp: num_major_snp,
                num_minor_snp: num_minor_snp,
                num_insertions: num_insertions,
                num_deletions: num_deletions,
                breadth_coverage: breadth_cov,
                depth_coverage: depth_cov,
                num_perfect_kmers: out_perfect_kmers,
                num_variant_kmers: out_variant_kmers,
                num_unmapped_kmers: out_unmapped_kmers,
            });
            variant_info.push((r1.to_string(), variants));
            
            
            cleanup_kmc_files(&args.output, &args.keep_kmer_counts);
            
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

fn cleanup_kmc_files(output_dir: &str, keep_counts: &bool) {
    if let Ok(entries) = fs::read_dir(output_dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            if let Some(name) = path.file_name().and_then(|n| n.to_str()) {
                if name.ends_with("kmc_pre")
                    || name.ends_with("kmc_suf")
                {
                    if let Err(e) = fs::remove_file(&path) {
                        warn!("Failed to delete {:?}: {}", path, e);
                    }
                }

                if !keep_counts {
                    if name.ends_with("_counts.txt") {
                        if let Err(e) = fs::remove_file(&path) {
                            warn!("Failed to delete {:?}: {}", path, e);
                        }
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
    let mut best_score: f64 = f64::NEG_INFINITY;

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

    if best_score == 0.0 {
        warn!("No genome had any unique kmers, suggesting sequencing depth below --min-kmers or some other issue")
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
    let mut best_score: f64 = f64::NEG_INFINITY;

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

    if best_score == 0.0 {
        warn!("No genome had any unique kmers, suggesting sequencing depth is below --min-kmers or some other issue")
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
            // SNP-only alignment: store the single ASCII base of each major substitution (indels skipped)
            if variant.af >= 0.5 && variant.ref_allele.len() == 1 && variant.alt_allele.len() == 1 {
                all_positions.insert((variant.seq.clone(), variant.pos), variant.ref_allele[0]);
                if let Some(sample_map) = sample_positions.get_mut(&sample.clone()) {
                    sample_map.insert((variant.seq.clone(), variant.pos), variant.alt_allele[0]);
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
            .map(|&b| b as char)
            .unwrap_or('N');
        ref_seq.push(ref_base);
    }
    writeln!(writer, ">{}", file_meta.name).unwrap();
    writeln!(writer, "{}", ref_seq).unwrap();

    for (sample_name, sample_map) in &sample_positions {
        let mut seq = String::with_capacity(positions.len());

        for pos_key in &positions {
            // If the sample has a variant at this position, use the alt base
            if let Some(&alt_base) = sample_map.get(pos_key) {
                seq.push(alt_base as char); // ASCII byte → char
            } else {
                // Otherwise fall back to reference base from all_positions
                let ref_base = all_positions
                    .get(pos_key)
                    .map(|&b| b as char)
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
    reconstructed: bool, // pileup is against the indel-corrected consensus (sets the reference_source column)
){
    info!("Writing output to pileup");

    let file_stem = clean_sample_id(read_output);

    let tsv_pileup = File::create(format!("{}/{}.tsv", args.output, file_stem)).unwrap_or_else(|e| {
        error!("{} | Failed to create tsv pileup file", e);
        std::process::exit(1);
    });
    let mut writer = BufWriter::new(tsv_pileup);

    writeln!(writer, "reference\treference_source\tindex\tref\tA\tC\tG\tT\ta\tc\tg\tt").unwrap();


    let file_meta = &viral_metadata.files[*best_genome_index as usize];

    // whether the pileup reflects the original reference or the indel-corrected consensus
    let reference_source = if reconstructed { "reconstructed" } else { "original" };

    for seq_entry in &file_meta.sequences {
        let seq = &seq_entry.name;

        let fwd = output.get(seq).expect("Could not match seq to fwd counts");
        let rev = output_rev.get(seq).expect("Could not match seq to rev counts");

        for (i, ref_base) in fwd.ref_bases.iter().enumerate() {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                seq,
                reference_source,
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

pub fn print_consensus(
    read_output: &String,
    args: &CallArgs,
    output: &DashMap<String, OutputData>,
    output_rev: &DashMap<String, OutputData>,
    viral_metadata: &ViralMetadata,
    best_genome_index: &u16,
    reconstructed_indels: &[ReconstructedIndel],
    major_variants: &[VCFRecord], // first-pass called variants; the major SNVs are applied to the consensus
){
    info!("Writing output to pileup");

    let file_stem = clean_sample_id(read_output);

    let fa_consensus = File::create(format!("{}/{}.fa", args.output, file_stem)).unwrap_or_else(|e| {
        error!("{} | Failed to create consensus fasta file", e);
        std::process::exit(1);
    });
    let mut writer = BufWriter::new(fa_consensus);

    let file_meta = &viral_metadata.files[*best_genome_index as usize];

    for seq_entry in &file_meta.sequences {
        let seq = &seq_entry.name;

        let fwd = output.get(seq).expect("Could not match seq to fwd counts");
        let rev = output_rev.get(seq).expect("Could not match seq to rev counts");

        // start from the reference: keep the reference base where covered, N where there is no coverage
        let mut consensus: Vec<u8> = Vec::with_capacity(fwd.ref_bases.len());
        for (i, &ref_base) in fwd.ref_bases.iter().enumerate() {
            let total_depth: u64 = (0..4).map(|b| fwd.counts[i][b] + rev.counts[i][b]).sum();
            consensus.push(if total_depth == 0u64 { b'N' } else { ref_base });
        }

        // apply the first-pass called major variants (AF >= 0.5 substitutions) for this sequence, so the
        // consensus only deviates from the reference at confident, filter-passing calls
        for v in major_variants {
            if v.seq == *seq && v.af >= 0.5 && v.ref_allele.len() == 1 && v.alt_allele.len() == 1 {
                if v.pos >= 1 && v.pos <= consensus.len() {
                    consensus[v.pos - 1] = v.alt_allele[0];
                }
            }
        }

        // splice the reconstructed indels for this sequence into the consensus (footprints share the
        // reference coordinate system, so they line up with the consensus buffer)
        let final_seq = splice_indels(&consensus, seq, reconstructed_indels);

        writeln!(writer, ">{}", seq).unwrap();
        for chunk in final_seq.chunks(DEFAULT_LINE_WRAP) {
            writer.write_all(chunk).unwrap();
            writer.write_all(b"\n").unwrap();
        }
        writer.write_all(b"\n").unwrap(); // blank line between sequences
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
        "filename\tselected_genome\tnum_major_snp\tnum_minor_snp\tnum_insertions\tnum_deletions\tbreadth_coverage\tdepth_coverage\tnum_perfect_kmers\tnum_variant_kmers\tnum_unmapped_kmers"
    ).unwrap();

    // rows
    for info in output_info {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{}\t{}\t{}",
            info.filename,
            info.selected_genome,
            info.num_major_snp,
            info.num_minor_snp,
            info.num_insertions,
            info.num_deletions,
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
    writeln!(writer, "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant type (SNP, MNV, INS, DEL)\">").unwrap();
    writeln!(writer, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO").unwrap();

    // sort records by contig (in the order contigs are declared above) then by position, so the
    // separately-appended indel records interleave with the SNVs instead of trailing at the end
    let contig_order: FxHashMap<&str, usize> = file_meta
        .sequences
        .iter()
        .enumerate()
        .map(|(i, s)| (s.name.as_str(), i))
        .collect();
    let mut sorted: Vec<&VCFRecord> = variants.iter().collect();
    sorted.sort_by(|a, b| {
        let ai = contig_order.get(a.seq.as_str()).copied().unwrap_or(usize::MAX);
        let bi = contig_order.get(b.seq.as_str()).copied().unwrap_or(usize::MAX);
        ai.cmp(&bi).then(a.pos.cmp(&b.pos))
    });

    for variant in sorted{
        let seq_out:&str = variant.seq.split_whitespace().next().unwrap_or("");
        writeln!(writer, "{}\t{}\t.\t{}\t{}\t.\tPASS\tDP={};AF={:.3};DP4={},{},{},{};SOR={:.3};TYPE={}",
            seq_out,
            variant.pos,
            String::from_utf8_lossy(&variant.ref_allele),
            String::from_utf8_lossy(&variant.alt_allele),
            variant.depth,
            variant.af,
            variant.fwd_ref, variant.rev_ref, variant.fwd_alt, variant.rev_alt,
            variant.sor,
            variant.var_type()
        ).unwrap()
    }
}

#[derive(Debug, Clone)]
pub struct VCFRecord{
    seq: String,
    pos: usize,
    ref_allele: Vec<u8>, // ASCII REF allele (length 1 for SNPs, multi-base for indels/MNVs)
    alt_allele: Vec<u8>, // ASCII ALT allele
    fwd_ref: u64,
    rev_ref: u64,
    fwd_alt: u64,
    rev_alt: u64,
    depth: u64,
    af: f64,
    sor: f64
}

impl VCFRecord {
    // VCF variant type inferred from allele lengths: equal-length-1 = SNP, equal-length>1 = MNV,
    // ALT longer = INS, REF longer = DEL
    fn var_type(&self) -> &'static str {
        match (self.ref_allele.len(), self.alt_allele.len()) {
            (r, a) if r == a && r == 1 => "SNP",
            (r, a) if r == a => "MNV",
            (r, a) if a > r => "INS",
            _ => "DEL",
        }
    }
}

// Build a VCF record for a reconstructed indel. REF is the original-reference footprint the allele
// replaced and ALT is the reconstructed allele. The two alleles drop out of opposite pileups, so each
// is measured where it maps cleanly: the REF allele (reads without the indel) from the first-pass
// pileup over the disrupted footprint, and the ALT allele (reads with the indel) from the second-pass
// (corrected consensus) pileup over the spliced-in allele. Each side takes the weakest (minimum total
// depth) position across its span, which lands in the disrupted region rather than the well-covered
// anchors. Those per-strand counts give DP4, AF, DP and the strand odds ratio. The shared flanks are
// then trimmed to a single anchor base and the start is reported in the original reference. Returns
// None if the indel can't be located/translated (e.g. footprint out of range or abutting another indel).
fn indel_to_vcf_record(
    indel: &ReconstructedIndel,
    orig_output: &DashMap<String, OutputData>,
    orig_output_rev: &DashMap<String, OutputData>,
    new_output: &DashMap<String, OutputData>,
    new_output_rev: &DashMap<String, OutputData>,
    coord_map: &CoordMap,
) -> Option<VCFRecord> {
    // REF allele: original reference bases over the 1-based inclusive footprint [ref_start, ref_end], plus
    // the first-pass reference-allele support. The reference allele is what still covers the disrupted
    // footprint in the first pass; the minimum-depth position lands in the dip (the deleted bases for a
    // deletion, or the junction dip for an insertion).
    let (ref_allele, ref_fwd, ref_rev) = {
        let orig = orig_output.get(&indel.seq)?;
        let orig_rev = orig_output_rev.get(&indel.seq)?;
        if indel.ref_start < 1 || indel.ref_end > orig.ref_bases.len() || indel.ref_start > indel.ref_end {
            return None;
        }
        let ref_allele = orig.ref_bases[indel.ref_start - 1..indel.ref_end].to_vec();

        let mut rep = indel.ref_start - 1;
        let mut min_total = u64::MAX;
        for p in (indel.ref_start - 1)..indel.ref_end {
            let total: u64 = (0..4).map(|b| orig.counts[p][b] + orig_rev.counts[p][b]).sum();
            if total < min_total {
                min_total = total;
                rep = p;
            }
        }
        let ref_base = nt_to_bits(orig.ref_bases[rep]) as usize;
        (ref_allele, orig.counts[rep][ref_base], orig_rev.counts[rep][ref_base])
    };

    // locate the allele in corrected coordinates: it begins just after the corrected position of the
    // reference base preceding the footprint (or at 0 when the footprint starts the sequence)
    let c_start = if indel.ref_start >= 2 {
        coord_map.orig_to_corrected.get(indel.ref_start - 2).copied().flatten()? + 1
    } else {
        0
    };

    // ALT allele: second-pass support over the spliced-in allele, again at the weakest spanning position
    let (alt_fwd, alt_rev) = {
        let fwd = new_output.get(&indel.seq)?;
        let rev = new_output_rev.get(&indel.seq)?;
        let c_end = (c_start + indel.allele.len()).min(fwd.counts.len());
        if c_start >= c_end {
            return None;
        }
        let mut rep = c_start;
        let mut min_total = u64::MAX;
        for p in c_start..c_end {
            let total: u64 = (0..4).map(|b| fwd.counts[p][b] + rev.counts[p][b]).sum();
            if total < min_total {
                min_total = total;
                rep = p;
            }
        }
        let alt_base = nt_to_bits(fwd.ref_bases[rep]) as usize;
        (fwd.counts[rep][alt_base], rev.counts[rep][alt_base])
    };

    let ref_total = ref_fwd + ref_rev;
    let alt_total = alt_fwd + alt_rev;
    let depth = ref_total + alt_total;
    let af = if depth > 0 { alt_total as f64 / depth as f64 } else { 0.0 };
    let sor = strand_odds_ratio(ref_fwd, ref_rev, alt_fwd, alt_rev);

    // trim the shared suffix entirely, then the shared prefix down to a single anchor base, shifting the
    // reported start past the trimmed prefix (standard compact indel representation)
    let mut ref_a = ref_allele;
    let mut alt_a = indel.allele.clone();

    let mut s = 0usize;
    while ref_a.len() - s > 1
        && alt_a.len() - s > 1
        && ref_a[ref_a.len() - 1 - s] == alt_a[alt_a.len() - 1 - s]
    {
        s += 1;
    }
    ref_a.truncate(ref_a.len() - s);
    alt_a.truncate(alt_a.len() - s);

    let mut p = 0usize;
    while p + 1 < ref_a.len() && p + 1 < alt_a.len() && ref_a[p] == alt_a[p] {
        p += 1;
    }
    let pos = indel.ref_start + p;
    let ref_a = ref_a[p..].to_vec();
    let alt_a = alt_a[p..].to_vec();

    if ref_a == alt_a {
        return None; // not an actual difference after trimming
    }

    Some(VCFRecord {
        seq: indel.seq.clone(),
        pos,
        ref_allele: ref_a,
        alt_allele: alt_a,
        fwd_ref: ref_fwd,
        rev_ref: ref_rev,
        fwd_alt: alt_fwd,
        rev_alt: alt_rev,
        depth,
        af,
        sor,
    })
}

#[derive(Debug, Clone)]
pub struct Noise{
    max: f64,
    mean: f64,
    std: f64,
}


pub fn get_baseline_noise(fwd: &OutputData, rev: &OutputData) -> Vec<Noise> {

    //window size and alppa, max table size for our streaming version of thompson tau
    let window_size = 100;
    let alpha = 0.001;
    let max_table_len = window_size / 10;

    let len = fwd.counts.len();
    
    // The output of what we want
    let mut baseline_noise = vec![Noise {max: 0.0, mean: 0.0, std: 0.0}; len];

    //storing the counts in the window (since max 3 alternative nucleotides the length is 3 * window length, which we can index based off of genome_pos % window_size)
    //we also keep in_max which is 0 or 1 depending on if it is in the max table
    let mut window_counts = vec![0.0; len*3];
    let mut in_max = vec![0; len*3];

    let mut maxes = vec![0.0; max_table_len];

    let mut n:usize = 0; //number of variant positions in the window (between 0-3*len)
    let mut s:f64 = 0.0; //sum of the variant positions
    let mut s2:f64 = 0.0; //sum^2 of variant positions
    let mut mu:f64; //rolling mean
    let mut var:f64; //rolling variance

    let half_window = window_size / 2;

    for i in 0..(len + half_window) {

        let base_pos = (i % window_size) * 3;
        
        
        let freqs: Vec<f64> = if i < len {
            // Combine counts across strands
            let mut counts: Vec<u64> = (0..4)
                .map(|b| fwd.counts[i][b] + rev.counts[i][b])
                .collect();
            counts.sort_unstable_by(|a, b| b.cmp(a));
            let total_depth: u64 = counts.iter().sum();
            if total_depth == 0 {
                vec![0.0; 4]
            } else {
                counts.iter().map(|&c| c as f64 / total_depth as f64).collect()
            }
        } else {
            vec![0.0; 4]
        };

        //loop through minor variants for that position and update n, s, s2
        for j in 1..4 {

            let idx = base_pos + (j - 1);

            //remove the old values
            let old = window_counts[idx];
            if old > 0.0 {
                n -= 1;
                s -= old;
                s2 -= old * old;

                //remove from max if there
                if in_max[idx] == 1 {
                    //remove from maxes table and shift others up
                    if let Some(pos) = maxes.iter().position(|&x| (x - old).abs() < 1e-12) {
                        for k in pos..(max_table_len - 1) {
                            maxes[k] = maxes[k+1];
                        }
                        maxes[max_table_len - 1] = 0.0;
                    }
                    in_max[idx] = 0;
                }
            }

            //update with new values
            let maf = freqs[j];
            if maf > 0.0 {
                n += 1;
                s += maf;
                s2 += maf * maf;

                //loop through max table from bottom and update
                for k in (0..max_table_len).rev() {
                    if maf > maxes[k] {
                        if k + 1 < max_table_len {
                            maxes[k + 1] = maxes[k];
                        }
                        maxes[k] = maf;
                    } else {
                        break;
                    }
                }
                in_max[idx] = 1;
            } else {
                in_max[idx] = 0;
                window_counts[idx] = 0.0;
            }

            window_counts[idx] = maf;

        }

        //update mu and var
        if n != 0 {
            mu = s / n as f64;
            var = (s2 / n as f64) - mu * mu;
        } else {
            mu = 0.0;
            var = 0.0;
        }

        //perform thompson tau for that window starting from largest in max:
        let mut curr_max_idx:usize = 0;
        let mut curr_n = n;
        let mut curr_s = s;
        let mut curr_s2 = s2;
        let mut curr_mu = mu;
        let mut curr_var = var;

        while curr_max_idx < max_table_len && maxes[curr_max_idx] != 0.0 {
            //calculate t from thompson tau test
            let candidate_outlier = maxes[curr_max_idx];
            let std = curr_var.sqrt();

            let tau = if curr_n > 2 {
                let df = (curr_n - 2) as f64;
                let t_crit = StudentsT::new(0.0, 1.0, df).unwrap()
                    .inverse_cdf(1.0 - alpha / (curr_n as f64));
                (t_crit * (curr_n as f64 - 1.0)) / ((curr_n as f64).sqrt() * ((curr_n as f64 - 2.0 + t_crit * t_crit).sqrt()))
            } else {
                f64::INFINITY
            };

            // see if current max - mu is greater than t - sqrt(var)
            // if so, adjust temp n, s, s2, then recalculate mu and var and go to next max. 
            // if not, set the baseline noise at i to curr_mu, curr_var, curr_max
            if (candidate_outlier as f64 - curr_mu).abs() > tau * std {
                curr_s -= candidate_outlier;
                curr_s2 -= candidate_outlier;
                curr_n -= 1;
                if curr_n > 0 {
                    curr_mu = curr_s as f64 / curr_n as f64;
                    curr_var = (curr_s2 as f64 / curr_n as f64) - curr_mu * curr_mu;
                } else {
                    curr_mu = 0.0;
                    curr_var = 0.0;
                }
                curr_max_idx += 1;
            } else {
                break;
            }

        }

        // update the noise levels
        if i >= half_window {
            let write_idx = i - half_window;
            if write_idx < len {
                baseline_noise[write_idx] = Noise{
                    max: maxes[curr_max_idx],
                    mean: curr_mu,
                    std: curr_var.sqrt()
                };
            }
        }
    }

    baseline_noise

}

// GATK strand odds ratio from a 2x2 strand table (ref/alt x fwd/rev), with +1 pseudocounts so it is
// always defined. Higher values indicate stronger strand bias.
fn strand_odds_ratio(ref_fwd: u64, ref_rev: u64, alt_fwd: u64, alt_rev: u64) -> f64 {
    let a = ref_fwd as f64 + 1.0; // ref fwd
    let b = ref_rev as f64 + 1.0; // ref rev
    let c = alt_fwd as f64 + 1.0; // alt fwd
    let d = alt_rev as f64 + 1.0; // alt rev

    let r = (a * d) / (b * c);
    let ref_ratio = a.min(b) / a.max(b);
    let alt_ratio = c.min(d) / c.max(d);

    (r + (1.0 / r)).ln() + ref_ratio.ln() - alt_ratio.ln()
}

// Recompute a first-pass major SNV against the second-pass (corrected consensus) pileup at its corrected
// position, so the reported depth/AF/strand counts reflect the improved mapping. REF stays the original
// reference base and ALT the major allele (now the consensus base); the variant is reported at its
// original-reference position. Returns None if the corrected position has no coverage (the caller then
// keeps the first-pass call).
fn major_snv_second_pass_record(
    seq: &str,
    orig_pos: usize,    // 1-based original-reference position
    corrected: usize,   // 0-based position in the corrected-consensus pileup
    ref_base: u8,       // ASCII original reference base
    alt_base: u8,       // ASCII major-variant allele (the consensus base)
    new_output: &DashMap<String, OutputData>,
    new_output_rev: &DashMap<String, OutputData>,
) -> Option<VCFRecord> {
    let fwd = new_output.get(seq)?;
    let rev = new_output_rev.get(seq)?;
    if corrected >= fwd.counts.len() {
        return None;
    }
    let total_depth: u64 = (0..4).map(|b| fwd.counts[corrected][b] + rev.counts[corrected][b]).sum();
    if total_depth == 0 {
        return None;
    }
    let rb = nt_to_bits(ref_base) as usize;
    let ab = nt_to_bits(alt_base) as usize;
    let fwd_ref = fwd.counts[corrected][rb];
    let rev_ref = rev.counts[corrected][rb];
    let fwd_alt = fwd.counts[corrected][ab];
    let rev_alt = rev.counts[corrected][ab];
    let af = (fwd_alt + rev_alt) as f64 / total_depth as f64;
    let sor = strand_odds_ratio(fwd_ref, rev_ref, fwd_alt, rev_alt);
    Some(VCFRecord {
        seq: seq.to_string(),
        pos: orig_pos,
        ref_allele: vec![ref_base],
        alt_allele: vec![alt_base],
        fwd_ref,
        rev_ref,
        fwd_alt,
        rev_alt,
        depth: total_depth,
        af,
        sor,
    })
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
    coord_map: Option<&FxHashMap<String, CoordMap>>, // second pass: translate corrected coords -> original reference
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

        let baseline_noise: Vec<Noise> = get_baseline_noise(&*fwd, &*rev);

        // for i in 100..120 {
        //     info!("{}", baseline_noise[i].mean);
        // }

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

            // calculate the depths, including those of fwd and reverse, find the min and max of the two for strand filtering purposes
            let row_total: Vec<u64> = (0..4)
                .map(|b| row[b] + row_rev[b])
                .collect();

            // tally coverage here, before any coord translation, so breadth/depth are computed the same
            // way regardless of pass -- over every covered position of whichever pileup table was passed
            // in. On the second pass this includes positions inside a reconstructed indel allele, which
            // are real covered positions of the corrected consensus.
            let total_depth = row_total.iter().sum();
            if total_depth == 0 {
                continue; //maybe replace down the line with logic to fix things
            } else {
                positions_covered += 1;
                total_coverage += total_depth as usize;
            }

            // on the second pass, `i` is a corrected-consensus coordinate; translate it back to the
            // original reference. Positions inside a reconstructed indel have no original coordinate, so
            // they are skipped for variant calling (the indel itself is reported as its own record).
            let pos = match coord_map {
                Some(maps) => match maps.get(seq).and_then(|m| m.corrected_to_orig.get(i).copied().flatten()) {
                    Some(o) => o + 1, // 1-based original coordinate
                    None => continue,
                },
                None => i + 1,
            };

            // loop through each base and variant call if not the reference base
            for alt_base in 0u8..4 {
                if alt_base == ref_base || row_total[alt_base as usize] == 0 { //skip if the reference base of if there isn't any depth at that position
                    continue;
                }

                // NEW Strand filter logic
                let mut sor = args.strand_odds_max + 1.0; //default is above the given so filtered if not able to calculate
                if strand_filter {
                    let a = row[ref_base as usize] as f64 + 1.0; //ref fwd
                    let b = row_rev[ref_base as usize] as f64 + 1.0; //ref rev
                    let c = row[alt_base as usize] as f64 + 1.0; //alt fwd
                    let d = row_rev[alt_base as usize] as f64 + 1.0; //alt rev

                    //calculate the difference between the strands
                    let ref_total = a + b + c + d;
                    let min_strand_depth = (a+c).min(b+d);
                    let min_strand_percent = min_strand_depth / ref_total;

                    //if the strand balance filter is on (default), then always do SOR filter
                    //if the strand balance filter is off, then only do SOR filter when 1 strand is > strand_balance_ratio of total depth (by default 10% of total depth, aka all cases where 1 strand is less than 10% of the total depth are ignored)
                    if (!args.no_strand_balance_filter) | ((args.no_strand_balance_filter) & (min_strand_percent >= args.strand_balance_ratio)) {

                        //Using GATK strand odds ratio
                        sor = strand_odds_ratio(row[ref_base as usize], row_rev[ref_base as usize], row[alt_base as usize], row_rev[alt_base as usize]);

                        // filter out if greater than strand odds ratio (default 8)
                        if sor > args.strand_odds_max {
                            continue;
                        }

                        // additional filtering for low kmer support across both strands
                        let c_k = count[alt_base as usize] as usize;
                        let d_k = count_rev[alt_base as usize] as usize;

                        if c_k < n_kmer_per_strand && d_k < n_kmer_per_strand {
                            continue;
                        }
                    } else {
                        sor = -1.0; //output SOR == -1.0 if the strand is unbalanced so not tested
                    }
                }

                //Get minor af (filter out if below reporting threshold)
                let alt_count = row_total[alt_base as usize];
                let af = alt_count as f64 / total_depth as f64;
                
                let y0: f64 = args.variant_multiplier;
                let p0: f64 = 0.5;
                let a: f64 = 0.03;
                let factor = y0 + p0 * a.powf(100.0 * af);

                if af < min_af || af < (factor.max(y0) * baseline_noise[i].max) {
                    continue;
                }


                // call major/minor variants
                if af >= 0.5 {
                    num_major_variants += 1;
                } else {
                    if total_depth < args.min_depth as u64 { // filter out minor variants when the total depth is too low
                        continue;
                    }
                    if alt_count < args.min_variant_depth as u64 { //filter out minor variants when variant depth is below threshold
                        continue;
                    }
                    num_minor_variants += 1;
                }

                results.push(VCFRecord {
                    seq: seq.clone(),
                    pos: pos,
                    ref_allele: vec![nucleotide_bits_to_char(ref_base as u64) as u8],
                    alt_allele: vec![nucleotide_bits_to_char(alt_base as u64) as u8],
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
    let mut depth_cov: f64 = 0.0;
    if positions_covered > 0 {
        depth_cov = total_coverage as f64 / positions_covered as f64;
    }
    info!("Sample breadth of coverage: {}, depth of coverage: {}", breadth_cov, depth_cov);
    info!("Called {} major variants, {} minor above maf = {}", num_major_variants, num_minor_variants, min_af);
    (results, num_major_variants, num_minor_variants, breadth_cov, depth_cov)

}

pub fn count_kmers_kmc(reads: &String, threads: &usize, args: &CallArgs) -> Result<(usize, usize, usize, usize), String> {
    let fastq_path = reads.clone();
    let file_stem = clean_sample_id(&fastq_path);


    let output_dir = Path::new(&args.output);
    let tmp_dir = output_dir.join(format!("tmp_{}", file_stem));
    trace!("{}, {}, {}", fastq_path, file_stem, tmp_dir.display());


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

    if !kmc_output.status.success() {
        let stderr = String::from_utf8_lossy(&kmc_output.stderr);
        return Err(format!("KMC3 counting failed: {}", stderr));
    }

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
    flanking_base_map: &DashMap<u64, [[u64; 4]; 2]>,
) -> FxHashMap<u16, (usize, usize, usize)> {
    let k = args.kmer;

    //to store the number of kmers that are mapped perfectly or with 1-edit distance after all of the chunks
    //key is the index of the file in ViralMetadata, values are number of perfectly mapped and 1-edit distance kmers
    let results: DashMap<u16, (usize, usize, usize)> = DashMap::new();
    
    let chunk_size = if ((kmers.len() / threads) as usize) < 10000 { max(kmers.len() / threads, 100) as usize } else { 10000 };

    let keep_mask: Option<Vec<bool>> = args.specify_pattern.as_ref().map(|sp| {
        parse_specify_pattern(sp, k).unwrap_or_else(|e| {
            error!("{}", e);
            std::process::exit(1);
        })
    });

    //let keep_mask = bucket_keep_mask(&args.bucket_pattern, args.bucket_stride);

    kmers.par_chunks(chunk_size).for_each(|chunk| {

        let mut local_counts: FxHashMap<u16, (usize, usize, usize)> = FxHashMap::default();

        // Thread-local output accumulators: (file_id, seq_id) -> position -> [base0..base3]
        // Avoids per-bucket DashMap lock; one DashMap write per sequence per chunk at merge time.
        type SeqKey = (u16, u8);
        let mut local_fwd_sum: FxHashMap<SeqKey, FxHashMap<usize, [u64; 4]>> = FxHashMap::default();
        let mut local_fwd_max: FxHashMap<SeqKey, FxHashMap<usize, [u64; 4]>> = FxHashMap::default();
        let mut local_rev_sum: FxHashMap<SeqKey, FxHashMap<usize, [u64; 4]>> = FxHashMap::default();
        let mut local_rev_max: FxHashMap<SeqKey, FxHashMap<usize, [u64; 4]>> = FxHashMap::default();

        let mut local_flanking: FxHashMap<u64, [[u64; 4]; 2]> = FxHashMap::default();

        for (kmer, n) in chunk {

            let (kmer_bin, rc) = canonical_kmer(kmer.as_bytes(), k);
            let buckets = assign_buckets(kmer_bin, k);

            //record the first and last bucket of the kmer with the base of the first/last character (canonical orientation), depth-weighted by n
            let rc_idx = rc as usize;
            let first_base = ((kmer_bin >> (2 * (k - 1))) & 0b11) as usize;
            let last_base = (kmer_bin & 0b11) as usize;
            local_flanking.entry(buckets[0]).or_insert([[0u64; 4]; 2])[rc_idx][first_base] += *n;
            local_flanking.entry(buckets[k - 1]).or_insert([[0u64; 4]; 2])[rc_idx][last_base] += *n;

            let filtered_buckets: Vec<u64> = if let Some(mask) = &keep_mask {
                buckets.iter().enumerate()
                    .filter(|(i, _)| mask[*i])
                    .map(|(_, &b)| b)
                    .collect()
            } else {
                let (start, end) = if args.use_full_kmer {
                    (0, buckets.len())
                } else {
                    let len = buckets.len();
                    if args.n_fixed * 2 + 1 >= len {
                        (0, 0)
                    } else {
                        (args.n_fixed, len - args.n_fixed - 1)
                    }
                };
                buckets[start..end].iter().enumerate()
                    .map(|(_, &b)| b)
                    .collect()
            };
            

            let num_buckets_perfect = filtered_buckets.len();
            let mut per_genome_bucket_hits: FxHashMap<u16, usize> = FxHashMap::default(); //number of hits to each genome (where the key is the index of the file) per kmer

            for &bucket in &filtered_buckets {

                if let Some(bucket_infos) = index.get(&bucket) {

                    for info in bucket_infos {
                        *per_genome_bucket_hits
                            .entry(info.file_id.clone())
                            .or_insert(0) += 1;

                        let genome_pos = info.location as usize;
                        let nuc_x = info.idx as usize;
                        let idx = genome_pos + nuc_x;
                        let key: SeqKey = (info.file_id, info.seq_id);

                        if info.canonical {
                            // pos = k - nuc_x - 1, so k - pos - 1 = nuc_x
                            let bit_idx = (((kmer_bin >> (2 * nuc_x)) & 0b11) ^ 0b11) as usize;
                            if rc {
                                local_fwd_sum.entry(key).or_default().entry(idx).or_insert([0u64; 4])[bit_idx] += 1;
                                let e = local_fwd_max.entry(key).or_default().entry(idx).or_insert([0u64; 4]);
                                if e[bit_idx] < *n { e[bit_idx] = *n; }
                            } else {
                                local_rev_sum.entry(key).or_default().entry(idx).or_insert([0u64; 4])[bit_idx] += 1;
                                let e = local_rev_max.entry(key).or_default().entry(idx).or_insert([0u64; 4]);
                                if e[bit_idx] < *n { e[bit_idx] = *n; }
                            }
                        } else {
                            // pos = nuc_x, so k - pos - 1 = k - nuc_x - 1
                            let bit_idx = ((kmer_bin >> (2 * (k - nuc_x - 1))) & 0b11) as usize;
                            if rc {
                                local_rev_sum.entry(key).or_default().entry(idx).or_insert([0u64; 4])[bit_idx] += 1;
                                let e = local_rev_max.entry(key).or_default().entry(idx).or_insert([0u64; 4]);
                                if e[bit_idx] < *n { e[bit_idx] = *n; }
                            } else {
                                local_fwd_sum.entry(key).or_default().entry(idx).or_insert([0u64; 4])[bit_idx] += 1;
                                let e = local_fwd_max.entry(key).or_default().entry(idx).or_insert([0u64; 4]);
                                if e[bit_idx] < *n { e[bit_idx] = *n; }
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
        // merge thread-local output buffers into shared DashMaps.
        
        for ((file_id, seq_id), pos_map) in local_fwd_sum {
            if let Some(maps) = output_maps.get(&file_id) {
                let seq = &viral_metadata.files[file_id as usize].sequences[seq_id as usize].name;
                if let Some(mut rec) = maps.2.get_mut(seq) {
                    for (pos, vals) in pos_map { for b in 0..4 { rec.counts[pos][b] += vals[b]; } }
                }
            }
        }
        for ((file_id, seq_id), pos_map) in local_fwd_max {
            if let Some(maps) = output_maps.get(&file_id) {
                let seq = &viral_metadata.files[file_id as usize].sequences[seq_id as usize].name;
                if let Some(mut rec) = maps.0.get_mut(seq) {
                    for (pos, vals) in pos_map { for b in 0..4 { if rec.counts[pos][b] < vals[b] { rec.counts[pos][b] = vals[b]; } } }
                }
            }
        }
        for ((file_id, seq_id), pos_map) in local_rev_sum {
            if let Some(maps) = output_maps.get(&file_id) {
                let seq = &viral_metadata.files[file_id as usize].sequences[seq_id as usize].name;
                if let Some(mut rec) = maps.3.get_mut(seq) {
                    for (pos, vals) in pos_map { for b in 0..4 { rec.counts[pos][b] += vals[b]; } }
                }
            }
        }
        for ((file_id, seq_id), pos_map) in local_rev_max {
            if let Some(maps) = output_maps.get(&file_id) {
                let seq = &viral_metadata.files[file_id as usize].sequences[seq_id as usize].name;
                if let Some(mut rec) = maps.1.get_mut(seq) {
                    for (pos, vals) in pos_map { for b in 0..4 { if rec.counts[pos][b] < vals[b] { rec.counts[pos][b] = vals[b]; } } }
                }
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

        // merge local flanking base counts into the global DashMap
        for (bucket, counts) in local_flanking {
            flanking_base_map
                .entry(bucket)
                .and_modify(|e| {
                    for r in 0..2 {
                        for b in 0..4 {
                            e[r][b] += counts[r][b];
                        }
                    }
                })
                .or_insert(counts);
        }
    });

    //ensure that each genome has 0s filled if nothing is mapped, otherwise error thrown downstream
    for (i, _) in viral_metadata.files.iter().enumerate() {
        let id = i as u16;
        results.entry(id).or_insert((0,0,0));
    }

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
