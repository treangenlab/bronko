use crate::cli::*;
use crate::io::*;
use crate::util::*;

use::log::*;

fn check_args(args: &CountArgs) {
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

    trace!("Test Trace");
    debug!("Test Debug");
    info!("Test Info");

    //Check kmer size to make sure it is odd and greater than 3
    if args.kmer % 2 != 1 || args.kmer < 3 {
        error!("Invalid kmer size, must be odd and >= 3");
        std::process::exit(1)
    }

    //check to see if input is fastq file
    if !check_fastq(&args.input) {
        error!("Input is not a fastq file (must be .fq(.gz)/.fastq(.gz)/.fnq(.gz))");
        std::process::exit(1)
    }
}

pub fn count(args: CountArgs) {
    // check to make sure none of the
    check_args(&args);

    match read_fastq(&args.input) {
        Ok((n_reads, n_bases)) => {
            info!("There are {} reads in your file", n_reads);
            info!("There are {} bases in your file.", n_bases);
        }
        Err(e) => match e {
            FastxError::IoError(io_err) => {
                error!("{}", io_err);
                std::process::exit(1);
            }
            FastxError::ParseError(parse_err) => {
                error!("Parse Error: {}", parse_err);
                std::process::exit(1);
            }
        },
    }
}
