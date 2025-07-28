use crate::cli::*;
use crate::util::*;

use log::*;
// use needletail::{FastxReader, parse_fastx_file};

// use std::collections::HashMap;

///
/// Struct for an entry in the hashmap, contains the kmer
/// 
// struct Kmer {
//     kmer: BitNuclKmer,
//     location: usize,
// }

fn check_args(args: &BuildArgs) {
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

    // trace!("Test Trace");
    // debug!("Test Debug");
    // info!("Test Info");

    //Check kmer size to make sure it is odd and greater than 3
    if args.kmer < 3 {
        error!("Invalid kmer size, must >= 3");
        std::process::exit(1)
    }

    //check to see if input is fastq file
    if !check_fasta(&args.genome) {
        error!("Input is not a fasta file (must be .fq(.gz)/.fastq(.gz)/.fnq(.gz))");
        std::process::exit(1)
    }

}


pub fn build(args: BuildArgs) {
    //check the arguments
    check_args(&args);

    info!("Starting build")

    // info!("Running Building");
    // let mut reader = parse_fastx_file(&args.genome).unwrap_or_else(|e| {
    //     error!("Fasta unable to be read: {}", e);
    //     std::process::exit(1);
    // });
    
    // while let Some(record) = reader.next() {
    //     let seqrec = record.unwrap_or_else(|e| {
    //         error!("Invalid record in {}: {}", &args.genome, e);
    //         std::process::exit(1);
    //     });
    // };    
    //     for i in 0..=seq.len().saturating_sub(k) {
    //         let kmer = &seq[i..i + k];

    //         // Get canonical k-mer and its binary representation
    //         if let Some(bit_kmer) = BitNuclKmer::from_text(kmer) {
    //             let canonical = bit_kmer.canonical();
    //             let binary = canonical.to_u64();

    //             // Store in HashSet (unique k-mers)
    //             unique_kmers.insert(binary);

    //             // Store in HashMap (k-mer counts)
    //             *kmer_counts.entry(binary).or_insert(0) += 1;
    //         } else {
    //             error!("Invalid k-mer at position {}", i);
    //         }
    //     }

    // }
    

    // Some nice statistics to have
    // let mut n_contigs: usize = 0;
    // let mut n_bases: usize = 0;
    // let mut n_kmers: usize = 0;

    // let mid: usize = &args.kmer / 2; //for finding the split kmer

    // let mut skmers: HashMap<Vec<u8>, [u32; 4]> = HashMap::new();

    // trace!("Building the skmer hashmap");
    // while let Some(record) = reader.next() {
    //     let seqrec = record.unwrap_or_else(|e| {
    //         error!("Invalid record in {}: {}", &args.genome, e);
    //         std::process::exit(1);
    //     });

    //     let norm_seq = seqrec.normalize(false);
    //     let seq = norm_seq.sequence();

    //     //for now just count reads and bases and return them
    //     n_contigs += 1;
    //     n_bases += seqrec.num_bases();
    // }

    // info!("Number of contigs in '{}': {}", &args.genome, n_contigs);
    // info!("Total bases in '{}': {}", &args.genome, n_bases);
    // info!("Kmer profile of genome ({} total, {} unique)", n_kmers, skmers.len());

    // let mut reader = parse_fastx_file(&args.input).unwrap_or_else(|e| {
    //     error!("Fastq unable to be read: {}", e);
    //     std::process::exit(1);
    // });

    // //FASTQ PARSING

    // let mut n_reads_fq = 0;
    // let mut n_bases_fq = 0;
    // let mut n_matching_kmers = 0;
    // let mut n_divergent_skmers = 0;

    // let mut skmers_divergent: HashMap<Vec<u8>, [u32; 4]> = HashMap::new();

    // while let Some(record) = reader.next() {
    //     let seqrec = record.unwrap_or_else(|e| {
    //         error!("Invalid record in {}: {}", &args.input, e);
    //         std::process::exit(1);
    //     });

    //     let norm_seq = seqrec.normalize(false);
    //     let seq = norm_seq.sequence();

    //     for kmer in seq.windows(args.kmer) {
    //         let (skmer, mid_base) = extract_kmer_parts(&kmer, mid);
    //         if let Some(entry) = skmers.get_mut(&skmer) {
    //             // skmer exists in skmers, update the corresponding nucleotide count
    //             match mid_base {
    //                 b'A' => entry[0] += 1,
    //                 b'C' => entry[1] += 1,
    //                 b'G' => entry[2] += 1,
    //                 b'T' => entry[3] += 1,
    //                 _ => (),
    //             }
    //             n_matching_kmers += 1;
    //         } else {
    //             let entry = skmers_divergent.entry(skmer).or_insert([0; 4]);
    //             match mid_base {
    //                 b'A' => entry[0] += 1,
    //                 b'C' => entry[1] += 1,
    //                 b'G' => entry[2] += 1,
    //                 b'T' => entry[3] += 1,
    //                 _ => (),
    //             }
    //             n_divergent_skmers += 1;
                
    //         }
    //     }
    //     n_reads_fq += 1;
    //     n_bases_fq += seqrec.num_bases();
    // }

    // info!("Number of reads in '{}': {}", &args.input, n_reads_fq);
    // info!("Total bases in '{}': {}", &args.input, n_bases_fq);
    // info!("Skmer profile of sequencing data ({} total, {} matched, {} divergent)", n_matching_kmers + n_divergent_skmers, n_matching_kmers, n_divergent_skmers);

    // // Loop through the skmer profile and print how many have divergent nucleotides
    // let threshold: u32 = 5; // TODO: will add as input parameter
    // let mut count_no_values = 0;
    // let mut count_one_value = 0;
    // let mut filtered_skmers: Vec<(&Vec<u8>, &[u32; 4])> = Vec::new();

    // for (skmer, counts) in &skmers {
    //     let above_threshold = counts.iter().filter(|&&count| count > threshold).count();

    //     if above_threshold > 1 {
    //         filtered_skmers.push((skmer, counts));
    //     } else if counts.iter().all(|&count| count == 0) {
    //         count_no_values += 1;
    //     } else if above_threshold == 1 {
    //         count_one_value += 1;
    //     }
    // }

    // info!("Skmer entries not found: {}/{}", count_no_values, n_kmers);
    // info!("Skmer entries with no mutational profile: {}", count_one_value);

    // info!("Potential intrahost variation sites: {}", filtered_skmers.len())
    

}

// Return the counts as a tuple}
