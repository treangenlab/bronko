use crate::cli::*;
use crate::util::*;
use crate::consts::*;
use crate::lcb::*;

use anyhow::{Result};
use anyhow::Error;

use log::*;

use rustc_hash::{FxHashMap};
use needletail::{parse_fastx_file};


use std::path::Path;

use rayon::prelude::*;


// Sequence metadata including name and length (individual segments/chromosomes)
pub struct SeqMeta {
    pub name: String,   //fasta header
    pub len: usize,     //length of sequence
    pub seq: Vec<u8>,     //sequence itself
}

/// Per File metadata with name and sequences (one genome)
pub struct FileMeta {
    pub name: String, //name of the file
    pub sequences: Vec<SeqMeta>, //sequences within the file (for segmented viruses, incomplete assemblies, etc)
}

/// All metadata (collection of all genomes in the index)
pub struct ViralMetadata {  
    pub files: Vec<FileMeta>,   //list of all files
    pub k: usize, //kmer size used in the index
}

#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct BucketInfo {
    pub file_id: u16, //index storing the filename 
    pub seq_id: u8, //index storing the sequence
    pub location: u32, //location of kmer in the sequence
    pub idx: u8, //location of bucket within kmer
    pub canonical: bool, //is the bucket from a canonical kmer
}

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


    //Check kmer size to make sure it is odd and greater than 3
    if args.kmer % 2 != 1 || args.kmer > MAX_KMER_SIZE || args.kmer < MIN_KMER_SIZE {
        error!("Invalid kmer size, must be odd and between [{}-{}]", MIN_KMER_SIZE, MAX_KMER_SIZE);
        std::process::exit(1)
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

}

pub fn build(args: BuildArgs) {

    //Check arguments for building index
    check_args(&args);

    //build the indexes
    let (ref_index, viral_metadata): (FxHashMap<u64, Vec<BucketInfo>>, ViralMetadata) = build_indexes(args.kmer, &args.genomes).unwrap_or_else(|e| {
        error!("{} | Reference failed to build", e);
        std::process::exit(1)
    });
    log_memory_usage(true, "Fasta files indexed successfully");
}

pub fn build_indexes(
    k: usize,
    genomes: &[String]
) -> Result<(FxHashMap<u64, Vec<BucketInfo>>, ViralMetadata), Error> {
    info!("Building indexes from fasta files");

    // Step 1: Each thread builds its own index + metadata
    let per_file: Vec<(FxHashMap<u64, Vec<BucketInfo>>, FileMeta)> = genomes
        .par_iter()
        .enumerate()
        .map(|(file_id, file_path)| {
            let mut reader = parse_fastx_file(file_path).unwrap_or_else(|e| {
                error!("{} | Failed to parse fasta file: {}", e, file_path);
                std::process::exit(1)
            });

            let file_name = Path::new(file_path)
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or("unknown")
                .to_string();

            let mut local_index: FxHashMap<u64, Vec<BucketInfo>> = FxHashMap::default();
            let mut sequences: Vec<SeqMeta> = Vec::new();

            let mut seq_id: u8 = 0;
            while let Some(record) = reader.next() {
                let record = record.unwrap_or_else(|e| {
                    error!("{} | Unable to read record in {}", e, file_path);
                    std::process::exit(1)
                });

                let seq = record.seq();
                let seq_name = String::from_utf8_lossy(record.id())
                    .split_whitespace()
                    .next()
                    .unwrap_or_default()
                    .to_string();
                let seq_len = seq.len();

                sequences.push(SeqMeta {
                    name: seq_name.clone(),
                    len: seq_len,
                    seq: seq.to_vec(), 
                });

                for i in 0..=seq_len.saturating_sub(k) {
                    let kmer = &seq[i..i + k];
                    let (kmer_bin, canonical) = canonical_kmer(kmer, k);
                    let buckets = assign_buckets(kmer_bin, k);

                    for (j, bucket_id) in buckets.iter().enumerate() {
                        local_index.entry(*bucket_id).or_default().push(BucketInfo {
                            file_id: file_id as u16,
                            seq_id,
                            location: i as u32,
                            idx: j as u8,
                            canonical: canonical,
                        });
                    }
                }

                seq_id += 1;
            }

            let file_meta = FileMeta {
                name: file_name,
                sequences,
            };

            (local_index, file_meta)
        })
        .collect();

    // Step 2: Merge all local indexes + collect metadata
    let mut global_index: FxHashMap<u64, Vec<BucketInfo>> = FxHashMap::default();
    let mut files: Vec<FileMeta> = Vec::with_capacity(per_file.len());

    for (local_index, file_meta) in per_file {
        for (bucket_id, mut entries) in local_index {
            global_index.entry(bucket_id).or_default().append(&mut entries);
        }
        files.push(file_meta);
    }

    Ok((global_index, ViralMetadata { files, k}))
}

// pub fn build_indexes(args: &BuildArgs) -> Result<(FxHashMap<u64, Vec<BucketInfo>>, DashMap<String, usize>), Error> {

//     info!("Building indexes from fasta files");
//     let k = args.kmer;
//     let index: Arc<Mutex<FxHashMap<u64, Vec<BucketInfo>>>> = Arc::new(Mutex::new(FxHashMap::default()));
//     let seq_info: DashMap<String, usize> = DashMap::new();

//     ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();

//     args.genomes.par_iter().for_each(|file_path|{
//         let mut reader = parse_fastx_file(file_path).unwrap_or_else(|e|{
//             error!("{} | Failed to parse fasta file: {}", e, file_path);
//             std::process::exit(1)
//         });

//         let file_name = Path::new(file_path).file_stem().and_then(|s| s.to_str()).unwrap_or("unknown").to_string();

//         while let Some(record) = reader.next() {
//             let record = record.unwrap_or_else(|e|{
//                 error!("{} | Unable to read record in {}", e, file_path);
//                 std::process::exit(1)
//             });

//             let seq: std::borrow::Cow<'_, [u8]> = record.seq();
//             let seq_name = String::from_utf8_lossy(record.id()).to_string().split_whitespace().next().unwrap_or_default().to_string();
//             let seq_len = seq.len();
            
//             seq_info.entry(seq_name.clone()).or_insert(seq_len);

//             // trace!("Identified {}", seq_name);
//             trace!("{} kmers in {}", seq.len().saturating_sub(k), seq_name);
//             for i in 0..=seq.len().saturating_sub(k){
//                 let kmer: &[u8] = &seq[i..i+k];
//                 let (kmer_bin, canonical) = canonical_kmer(kmer, k);
//                 let buckets = assign_buckets(kmer_bin, k);

//                 for (j, bucket_id) in buckets.iter().enumerate(){
//                     let bucket = BucketInfo {
//                         file_name: file_name.clone(), //might want to change this in the future as it is kind of redundant
//                         seq_name: seq_name.clone(),
//                         location: i, //should be able to combine location + idx into single location
//                         idx: j,
//                         ref_base: nt_to_bits(kmer[j]) as u8, //this also is redundant
//                         canonical: canonical //this is maybe the only important thing aside from location
//                     };

//                     let mut locked = index.lock().unwrap();
//                     locked.entry(*bucket_id).or_default().push(bucket);
//                 }
//             }
//         }
//         info!("Processed sample")

//     });

//     Ok((
//         Arc::try_unwrap(index)
//             .unwrap()
//             .into_inner()
//             .unwrap(),
//         seq_info
//     ))
// }
