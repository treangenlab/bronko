pub fn check_fastq(file: &str) -> bool {
    if file.ends_with(".fq")
        || file.ends_with(".fastq")
        || file.ends_with(".fq.gz")
        || file.ends_with("fastq.gz")
        || file.ends_with("fnq")
        || file.ends_with("fnq.gz")
    {
        return true;
    }
    return false;
}

pub fn check_fasta(file: &str) -> bool {
    if file.ends_with(".fa")
        || file.ends_with(".fasta")
        || file.ends_with(".fa.gz")
        || file.ends_with("fasta.gz")
        || file.ends_with("fna")
        || file.ends_with("fna.gz")
    {
        return true;
    }
    return false;
}