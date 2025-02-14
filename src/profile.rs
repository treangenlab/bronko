use crate::cli::*;
use crate::util::*;

use log::*;
use needletail::{parse_fastx_file, Sequence};

use std::collections::HashMap;

fn check_args(args: &ProfileArgs) {
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
    if args.kmer % 2 != 1 || args.kmer < 3 {
        error!("Invalid kmer size, must be odd and >= 3");
        std::process::exit(1)
    }

    //check to see if input is fastq file
    if !check_fasta(&args.input) {
        error!("Input is not a fastq file (must be .fq(.gz)/.fastq(.gz)/.fnq(.gz))");
        std::process::exit(1)
    }
}

pub fn profile(args: ProfileArgs) {
    //check the arguments
    check_args(&args);

    info!("Running Profiling");
    let mut reader = parse_fastx_file(&args.input).unwrap_or_else(|e| {
        error!("Fasta unable to be read: {}", e);
        std::process::exit(1);
    });

    // Some nice statistics to have
    let mut n_reads: usize = 0;
    let mut n_bases: usize = 0;
    let mut n_kmers: usize = 0;

    let mid: usize = &args.kmer / 2; //for finding the split kmer

    let mut skmers: HashMap<Vec<u8>, [u32; 4]> = HashMap::new();

    trace!("Building the skmer hashmap");
    while let Some(record) = reader.next() {
        let seqrec = record.unwrap_or_else(|e| {
            error!("Invalid record in {}: {}", &args.input, e);
            std::process::exit(1);
        });

        let norm_seq = seqrec.normalize(false);
        let seq = norm_seq.sequence();

        for kmer in seq.windows(args.kmer) {
            let skmer: Vec<u8> = [&kmer[..mid], &kmer[mid + 1..]].concat();
            let mid = &kmer[mid];

            let entry = skmers.entry(skmer).or_insert([0; 4]);

            match mid {
                b'A' => entry[0] += 1,
                b'C' => entry[1] += 1,
                b'G' => entry[2] += 1,
                b'T' => entry[3] += 1,
                _ => (),
            }

            n_kmers += 1;
        }

        //for now just count reads and bases and return them
        n_reads += 1;
        n_bases += seqrec.num_bases();
    }

    info!("Number of contigs in '{}': {}", &args.input, n_reads);
    info!("Total bases in '{}': {}", &args.input, n_bases);
    info!("Kmer profile ({} total, {} unique)", n_kmers, skmers.len());
}

// Return the counts as a tuple}
