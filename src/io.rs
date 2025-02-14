use needletail::{errors, parse_fastx_file};

#[derive(Debug)]
pub enum FastxError {
    IoError(errors::ParseError),
    ParseError(String),
}

impl From<errors::ParseError> for FastxError {
    fn from(error: errors::ParseError) -> Self {
        FastxError::IoError(error)
    }
}

impl From<String> for FastxError {
    fn from(error: String) -> Self {
        FastxError::ParseError(error)
    }
}

pub fn read_fastq(file: &str) -> Result<(usize, usize), FastxError> {
    let mut n_bases: usize = 0;
    let mut n_reads: usize = 0;

    // Try to parse the FASTQ file and handle errors
    let mut reader = parse_fastx_file(file).map_err(|e| FastxError::IoError(e))?;

    while let Some(record) = reader.next() {
        let seqrec =
            record.map_err(|e| FastxError::ParseError(format!("Invalid record: {}", e)))?;

        //for now just count reads and bases and return them
        n_reads += 1;
        n_bases += seqrec.num_bases();
    }

    // Return the counts as a tuple
    Ok((n_reads, n_bases))
}
