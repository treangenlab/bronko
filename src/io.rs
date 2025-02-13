use needletail::{parse_fastx_file, errors};

#[derive(Debug)]
pub enum FastqError {
    IoError(errors::ParseError),
    ParseError(String)
}

impl From<errors::ParseError> for FastqError {
    fn from(error: errors::ParseError) -> Self {
        FastqError::IoError(error)
    }
}

impl From<String> for FastqError {
    fn from(error: String) -> Self {
        FastqError::ParseError(error)
    }
}

// Assuming `parse_fastx_file` and other functions are defined elsewhere
pub fn read_fastq(file: &str) -> Result<(usize, usize), FastqError> {
    let mut n_bases: usize = 0;
    let mut n_reads: usize = 0;

    // Try to parse the FASTQ file and handle errors
    let mut reader = parse_fastx_file(file).map_err(|e| FastqError::IoError(e))?;

    while let Some(record) = reader.next() {
        let seqrec = record.map_err(|e| FastqError::ParseError(format!("Invalid record: {}", e)))?;

        //for now just count reads and bases and return them
        n_reads += 1;
        n_bases += seqrec.num_bases();

    }

    // Return the counts as a tuple
    Ok((n_reads, n_bases))
}


