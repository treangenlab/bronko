use crate::cli::*;
use crate::util::*;
use crate::consts::*;
use log::*;

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
    
}


// pub fn build_indexes(args: &QueryArgs) -> Result<(FxHashMap<u64, Vec<BucketInfo>>, DashMap<String, usize>), Error> {

//     info!("Building indexes from fasta files");
//     let k = args.kmer;
//     let index: Arc<Mutex<FxHashMap<u64, Vec<BucketInfo>>>> = Arc::new(Mutex::new(FxHashMap::default()));
//     let seq_info: DashMap<String, usize> = DashMap::new();

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

//     });

//     Ok((
//         Arc::try_unwrap(index)
//             .unwrap()
//             .into_inner()
//             .unwrap(),
//         seq_info
//     ))
// }
