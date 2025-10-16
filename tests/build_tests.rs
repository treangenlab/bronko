use assert_cmd::Command;

#[test]
fn test_build_4_genomes() {
    let output_dir = tempfile::tempdir().unwrap();
    let output_path = output_dir.path().join("bronko");

    let mut cmd = Command::cargo_bin("bronko").unwrap();
    cmd.arg("build")
        .arg("-g")
        .arg("test_data/4_sarscov2/wuhan_ref.fasta")
        .arg("test_data/4_sarscov2/OM223929.1.fasta")
        .arg("test_data/4_sarscov2/ON765678.1.fasta")
        .arg("test_data/4_sarscov2/PX392231.1.fasta")
        .arg("-o")
        .arg(&output_path)
        .assert()
        .success()
        .code(0);
}

#[test]
fn test_build_single_genome_and_k19() {
    let output_dir = tempfile::tempdir().unwrap();
    let output_path = output_dir.path().join("hpv_index");

    Command::cargo_bin("bronko")
        .unwrap()
        .args([
            "build",
            "-g",
            "test_data/HPV16.fa",
            "-k",
            "19",
            "-o",
            output_path.to_str().unwrap(),
        ])
        .assert()
        .success()
        .code(0);
}

#[test]
fn test_build_single_genome() {

    Command::cargo_bin("bronko")
        .unwrap()
        .args([
            "build",
            "-g",
            "test_data/HPV16.fa",
        ])
        .assert()
        .success()
        .code(0);
}