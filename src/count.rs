use log::*;
use crate::cli::*;

fn check_args(args: &CountArgs){

    let output_level;
    if args.verbose {
        output_level = log::LevelFilter::Trace;
    } else if args.debug {
        output_level = log::LevelFilter::Debug;
    } else {
        output_level = log::LevelFilter::Info;
    }

    simple_logger::SimpleLogger::new().with_level(output_level).init().unwrap();

    trace!("Test Trace");
    debug!("Test Debug");
    info!("Test Info");

    //Check kmer size to make sure it is odd and greater than 3
    if args.kmer % 2 != 1 || args.kmer < 3 {
        error!("Invalid kmer size, must be odd and >= 3");
        std::process::exit(1)
    }


}

pub fn count(args: CountArgs){

    check_args(&args);
    
}