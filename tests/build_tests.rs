use assert_cmd::Command;

#[test]
fn build_tests() {
    let output_dir = tempfile::tempdir().unwrap();
    let output_path = output_dir.path().join("bronko");

    let mut cmd = Command::cargo_bin("bronko").unwrap();
    cmd.arg("build")
        .arg("-g")
        .arg("test_data/4_sarscov2/wuhan_ref.fasta")
        .arg("test_data/4_sarscov2/OM223929.1.fasta")
        .arg("test_data/4_sarscov2/ON765678.1.fasta")
        .arg("test_data/4_sarscov2/PX392231.1.fasta")
        .arg("-t")
        .arg("2")
        .arg("-o")
        .arg(&output_path)
        .assert()
        .success()
        .code(0);

    let mut cmd = Command::cargo_bin("bronko").unwrap();
    cmd.arg("build")
        .arg("-g")
        .arg("test_data/HPV16.fa")
        .arg("-k")
        .arg("19")
        .arg("-t")
        .arg("2")
        .arg("-o")
        .arg(&output_path)
        .assert()
        .success()
        .code(0);

    let mut cmd = Command::cargo_bin("bronko").unwrap();
    cmd.arg("build")
        .arg("-g")
        .arg("test_data/HPV16.fa")
        .arg("-t")
        .arg("2")
        .arg("-o")
        .arg(&output_path)
        .assert()
        .success()
        .code(0);

}

