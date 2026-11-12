//! End-to-end tests: run the program on the tiny made-up transcriptome in
//! `tests/data` and check its output files.

use std::collections::HashMap;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_rnaseq-readsim");

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("data")
        .join(name)
}

/// A scratch directory that is removed when the test ends.
struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("rnaseq-readsim-test-{}-{name}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        Scratch(dir)
    }

    fn path(&self, name: &str) -> String {
        self.0.join(name).to_string_lossy().into_owned()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn run(args: &[&str]) -> Output {
    Command::new(BIN)
        .args(args)
        .output()
        .expect("the program starts")
}

fn simulate(prefix: &str, extra: &[&str]) -> Output {
    let transcripts = data("transcripts.fa");
    let abundance = data("abundance.tsv");
    let mut args = vec![
        "simulate",
        "--transcripts",
        transcripts.to_str().unwrap(),
        "--abundance",
        abundance.to_str().unwrap(),
        "--output-prefix",
        prefix,
        "--quiet",
    ];
    args.extend_from_slice(extra);
    run(&args)
}

fn assert_success(out: &Output) {
    assert!(
        out.status.success(),
        "the program failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}

fn read_bytes(path: &str) -> Vec<u8> {
    std::fs::read(path).unwrap_or_else(|e| panic!("cannot read {path}: {e}"))
}

fn gunzip(path: &str) -> Vec<u8> {
    let mut out = Vec::new();
    flate2::read::MultiGzDecoder::new(&read_bytes(path)[..])
        .read_to_end(&mut out)
        .unwrap_or_else(|e| panic!("{path} is not valid gzip: {e}"));
    out
}

struct Record {
    name: String,
    comment: String,
    seq: String,
    qual: String,
}

fn fastq(text: &[u8]) -> Vec<Record> {
    let text = std::str::from_utf8(text).expect("FASTQ is ASCII");
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len() % 4, 0, "FASTQ has a multiple of four lines");
    lines
        .chunks(4)
        .map(|c| {
            assert!(c[0].starts_with('@'), "header: {}", c[0]);
            assert_eq!(c[2], "+");
            assert_eq!(c[1].len(), c[3].len());
            let (name, comment) = c[0][1..].split_once(' ').expect("header has a comment");
            Record {
                name: name.to_owned(),
                comment: comment.to_owned(),
                seq: c[1].to_owned(),
                qual: c[3].to_owned(),
            }
        })
        .collect()
}

fn transcripts() -> HashMap<String, String> {
    let text = std::fs::read_to_string(data("transcripts.fa")).unwrap();
    let mut map = HashMap::new();
    let mut id = String::new();
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            id = header.split('|').next().unwrap().to_owned();
            map.insert(id.clone(), String::new());
        } else {
            map.get_mut(&id).unwrap().push_str(line.trim());
        }
    }
    map
}

fn reverse_complement(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            'T' => 'A',
            _ => 'N',
        })
        .collect()
}

/// A profile in which every base has quality 93: no substitutions in
/// practice, so reads can be compared with the transcripts exactly.
const PERFECT_PROFILE: &str = r#"{
  "format": "rnaseq-readsim-quality-profile",
  "version": 1,
  "description": "every base Q93",
  "reads": [ { "qualities": [93], "cycles": [[1]] } ]
}"#;

#[test]
fn default_run_writes_consistent_files() {
    let scratch = Scratch::new("default");
    let prefix = scratch.path("sample");
    assert_success(&simulate(&prefix, &["-n", "5000", "--threads", "2"]));

    let r1 = fastq(&gunzip(&format!("{prefix}_1.fastq.gz")));
    let r2 = fastq(&gunzip(&format!("{prefix}_2.fastq.gz")));
    assert_eq!(r1.len(), 5000);
    assert_eq!(r2.len(), 5000);
    let mut names = std::collections::HashSet::new();
    for (a, b) in r1.iter().zip(&r2) {
        assert_eq!(a.name, b.name);
        assert_eq!(a.comment, "1:N:0:CCGCGGTT+CTAGCGCT");
        assert_eq!(b.comment, "2:N:0:CCGCGGTT+CTAGCGCT");
        assert!(
            names.insert(a.name.clone()),
            "duplicate read name {}",
            a.name
        );
        // @<instrument>:<run>:<flowcell>:<lane>:<tile>:<x>:<y>
        let fields: Vec<&str> = a.name.split(':').collect();
        assert_eq!(fields.len(), 7);
        assert_eq!(&fields[..4], &["A00001", "1", "H7SIMDSXX", "1"]);
        for f in &fields[4..] {
            f.parse::<u32>().expect("tile, x and y are numbers");
        }
        // the truth is not in the read names
        assert!(!a.name.contains("TOY"));
        for r in [a, b] {
            assert_eq!(r.seq.len(), 100);
            assert!(r.seq.bytes().all(|c| b"ACGTN".contains(&c)));
            assert!(r.qual.bytes().all(|c| b",:F".contains(&c)), "{}", r.qual);
        }
    }

    // truth table: one row per transcript, counts add up, zero stays zero
    let truth = std::fs::read_to_string(format!("{prefix}.truth.tsv")).unwrap();
    let mut lines = truth.lines();
    assert_eq!(
        lines.next().unwrap(),
        "transcript_id\tlength\teffective_length\tabundance\texpected_fragments\tfragments"
    );
    let rows: Vec<Vec<&str>> = lines.map(|l| l.split('\t').collect()).collect();
    assert_eq!(rows.len(), 8);
    let total: u64 = rows.iter().map(|r| r[5].parse::<u64>().unwrap()).sum();
    assert_eq!(total, 5000);
    let expected: f64 = rows.iter().map(|r| r[4].parse::<f64>().unwrap()).sum();
    assert!((expected - 5000.0).abs() < 0.01);
    let zero = rows.iter().find(|r| r[0] == "TOY0008.1").unwrap();
    assert_eq!((zero[3], zero[5]), ("0", "0"));
    assert_eq!(rows[0][0], "TOY0001.1");
    assert_eq!(rows[0][1], "2400");

    // report
    let report: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(format!("{prefix}.report.json")).unwrap())
            .unwrap();
    assert_eq!(report["tool"], "rnaseq-readsim");
    assert_eq!(report["version"], env!("CARGO_PKG_VERSION"));
    assert_eq!(report["parameters"]["seed"], 1);
    assert_eq!(report["parameters"]["library_type"], "reverse");
    assert_eq!(report["parameters"]["read_length"], 100);
    assert_eq!(report["parameters"]["quality_profile"], "built-in");
    assert_eq!(report["summary"]["fragments"], 5000);
    assert_eq!(report["input"]["transcripts"], 8);
    let short = report["summary"]["fragments_shorter_than_read"]
        .as_u64()
        .unwrap();
    assert!(short > 0 && short < 5000);
    assert!(report["summary"]["fragment_length_mean"].as_f64().unwrap() > 100.0);
    // the per-read truth is off by default
    assert!(!Path::new(&format!("{prefix}.reads.tsv.gz")).exists());

    // adapter read-through is visible in the reads
    let with_adapter = r1
        .iter()
        .filter(|r| r.seq.contains("AGATCGGAAGAGC"))
        .count() as u64;
    assert!(with_adapter > 0 && with_adapter <= short);
}

#[test]
fn output_does_not_depend_on_the_number_of_threads() {
    let scratch = Scratch::new("threads");
    // more than four chunks, so that threads really interleave
    let common = ["-n", "20000", "--per-read-truth", "--seed", "99"];
    let mut files: Vec<Vec<Vec<u8>>> = Vec::new();
    for threads in ["1", "2", "5"] {
        let prefix = scratch.path(&format!("t{threads}"));
        let mut args = common.to_vec();
        args.extend_from_slice(&["--threads", threads]);
        assert_success(&simulate(&prefix, &args));
        files.push(
            ["_1.fastq.gz", "_2.fastq.gz", ".reads.tsv.gz", ".truth.tsv"]
                .iter()
                .map(|ext| read_bytes(&format!("{prefix}{ext}")))
                .collect(),
        );
        // the report differs only in the file names it lists
        let report = std::fs::read_to_string(format!("{prefix}.report.json")).unwrap();
        files
            .last_mut()
            .unwrap()
            .push(report.replace(&prefix, "PREFIX").into_bytes());
    }
    assert_eq!(files[0], files[1], "1 and 2 threads differ");
    assert_eq!(files[0], files[2], "1 and 5 threads differ");

    // plain output: same reads as the gzip output, and again thread-independent
    let mut plain: Vec<Vec<u8>> = Vec::new();
    for threads in ["1", "3"] {
        let prefix = scratch.path(&format!("u{threads}"));
        let mut args = common.to_vec();
        args.extend_from_slice(&["--threads", threads, "--uncompressed"]);
        assert_success(&simulate(&prefix, &args));
        plain.push(read_bytes(&format!("{prefix}_1.fastq")));
        assert_eq!(
            read_bytes(&format!("{prefix}_2.fastq")),
            gunzip(&scratch.path("t1_2.fastq.gz"))
        );
        assert_eq!(
            read_bytes(&format!("{prefix}.reads.tsv")),
            gunzip(&scratch.path("t1.reads.tsv.gz"))
        );
    }
    assert_eq!(plain[0], plain[1]);
    assert_eq!(plain[0], gunzip(&scratch.path("t1_1.fastq.gz")));

    // another seed gives other reads; the same seed gives the same file again
    let other = scratch.path("other");
    assert_success(&simulate(
        &other,
        &["-n", "20000", "--seed", "100", "--uncompressed"],
    ));
    assert_ne!(read_bytes(&format!("{other}_1.fastq")), plain[0]);
    let again = scratch.path("again");
    assert_success(&simulate(
        &again,
        &["-n", "20000", "--seed", "99", "--uncompressed"],
    ));
    assert_eq!(read_bytes(&format!("{again}_1.fastq")), plain[0]);
}

#[test]
fn reads_agree_with_the_per_read_truth() {
    let scratch = Scratch::new("truth");
    let profile = scratch.path("perfect.json");
    std::fs::write(&profile, PERFECT_PROFILE).unwrap();
    let tx = transcripts();
    let adapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGCGGTTATCTCGTATGCCGTCTTCTGCTTG";
    let adapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTCTAGCGCTGTGTAGATCTCGGTGGTCGCCGTATCATT";

    for (library, salmon_like) in [("reverse", "-"), ("forward", "+"), ("unstranded", "")] {
        let prefix = scratch.path(library);
        assert_success(&simulate(
            &prefix,
            &[
                "-n",
                "4000",
                "--library-type",
                library,
                "--quality-profile",
                &profile,
                "--per-read-truth",
                "--adapter-fill",
                "polyA",
            ],
        ));
        let r1 = fastq(&gunzip(&format!("{prefix}_1.fastq.gz")));
        let r2 = fastq(&gunzip(&format!("{prefix}_2.fastq.gz")));
        let truth = String::from_utf8(gunzip(&format!("{prefix}.reads.tsv.gz"))).unwrap();
        let mut lines = truth.lines();
        assert_eq!(
            lines.next().unwrap(),
            "read_name\ttranscript_id\tstart\tfragment_length\tread1_strand"
        );
        let rows: Vec<Vec<&str>> = lines.map(|l| l.split('\t').collect()).collect();
        assert_eq!(rows.len(), 4000);
        let mut per_transcript: HashMap<String, u64> = HashMap::new();
        let mut plus = 0;
        let mut read_through = 0;
        let mut past_adapter = 0;
        for ((a, b), row) in r1.iter().zip(&r2).zip(&rows) {
            assert_eq!(a.name, row[0]);
            assert!(a.qual.bytes().all(|c| c == b'~'));
            let transcript = &tx[row[1]];
            let start: usize = row[2].parse::<usize>().unwrap() - 1;
            let len: usize = row[3].parse().unwrap();
            let fragment = &transcript[start..start + len];
            let rc = reverse_complement(fragment);
            if !salmon_like.is_empty() {
                assert_eq!(row[4], salmon_like);
            }
            plus += (row[4] == "+") as u32;
            let (insert1, insert2) = if row[4] == "+" {
                (fragment, rc.as_str())
            } else {
                (rc.as_str(), fragment)
            };
            for (read, insert, adapter) in
                [(&a.seq, insert1, adapter1), (&b.seq, insert2, adapter2)]
            {
                let n = len.min(100);
                assert_eq!(&read[..n], &insert[..n]);
                let k = adapter.len().min(100 - n);
                assert_eq!(&read[n..n + k], &adapter[..k]);
                assert!(read[n + k..].bytes().all(|c| c == b'A'));
            }
            read_through += (len < 100) as u32;
            past_adapter += (len + adapter1.len() < 100) as u32;
            *per_transcript.entry(row[1].to_owned()).or_default() += 1;
        }
        assert!(
            read_through > 100,
            "{read_through} fragments shorter than the read"
        );
        // no fragment is shorter than 40 bases, and 40 + 66 adapter bases
        // exceed the read length: no read reaches the filler
        assert_eq!(past_adapter, 0);
        if library == "unstranded" {
            assert!(
                (1800..2200).contains(&plus),
                "{plus} of 4000 on the + strand"
            );
        }

        // the per-transcript truth table agrees with the per-read table
        let table = std::fs::read_to_string(format!("{prefix}.truth.tsv")).unwrap();
        for line in table.lines().skip(1) {
            let f: Vec<&str> = line.split('\t').collect();
            let n: u64 = f[5].parse().unwrap();
            assert_eq!(
                n,
                per_transcript.get(f[0]).copied().unwrap_or(0),
                "{}",
                f[0]
            );
        }
    }
}

#[test]
fn reads_run_past_a_short_adapter_into_the_filler() {
    let scratch = Scratch::new("filler");
    let profile = scratch.path("perfect.json");
    std::fs::write(&profile, PERFECT_PROFILE).unwrap();
    let prefix = scratch.path("s");
    assert_success(&simulate(
        &prefix,
        &[
            "-n",
            "2000",
            "--quality-profile",
            &profile,
            "--adapter1",
            "ACGT{i7}",
            "--adapter2",
            "TTT{i5}",
            "--index1",
            "CC",
            "--index2",
            "",
            "--per-read-truth",
            "--uncompressed",
        ],
    ));
    let r1 = fastq(&read_bytes(&format!("{prefix}_1.fastq")));
    let r2 = fastq(&read_bytes(&format!("{prefix}_2.fastq")));
    let truth = std::fs::read_to_string(format!("{prefix}.reads.tsv")).unwrap();
    let mut seen = 0;
    for ((a, b), line) in r1.iter().zip(&r2).zip(truth.lines().skip(1)) {
        assert_eq!(a.comment, "1:N:0:CC");
        let len: usize = line.split('\t').nth(3).unwrap().parse().unwrap();
        if len + 6 < 100 {
            seen += 1;
            assert_eq!(&a.seq[len..len + 6], "ACGTCC");
            assert!(a.seq[len + 6..].bytes().all(|c| c == b'G'));
            assert_eq!(&b.seq[len..len + 3], "TTT");
            assert!(b.seq[len + 3..].bytes().all(|c| c == b'G'));
        }
    }
    assert!(seen > 50, "only {seen} reads ran past the adapter");
    let report: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(format!("{prefix}.report.json")).unwrap())
            .unwrap();
    assert_eq!(report["summary"]["fragments_with_read1_past_adapter"], seen);
    assert_eq!(report["resolved_adapter1"], "ACGTCC");
    assert_eq!(report["resolved_adapter2"], "TTT");
}

#[test]
fn a_learned_profile_reproduces_the_qualities_it_was_learned_from() {
    let scratch = Scratch::new("profile");
    let first = scratch.path("first");
    assert_success(&simulate(&first, &["-n", "30000", "--read-length", "60"]));
    let profile = scratch.path("learned.json");
    let out = run(&[
        "profile",
        "-1",
        &format!("{first}_1.fastq.gz"),
        "-2",
        &format!("{first}_2.fastq.gz"),
        "-o",
        &profile,
        "--quiet",
    ]);
    assert_success(&out);
    let json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&profile).unwrap()).unwrap();
    assert_eq!(json["format"], "rnaseq-readsim-quality-profile");
    assert_eq!(json["reads"].as_array().unwrap().len(), 2);
    assert_eq!(
        json["reads"][0]["qualities"],
        serde_json::json!([11, 25, 37])
    );
    assert_eq!(json["reads"][0]["cycles"].as_array().unwrap().len(), 60);
    assert_eq!(
        json["reads"][0]["transitions"].as_array().unwrap().len(),
        59
    );

    let second = scratch.path("second");
    assert_success(&simulate(
        &second,
        &[
            "-n",
            "30000",
            "--read-length",
            "60",
            "--quality-profile",
            &profile,
            "--seed",
            "5",
        ],
    ));
    // fraction of Q37 and of Q11-after-Q11 in both runs
    let summarise = |prefix: &str| -> (f64, f64) {
        let reads = fastq(&gunzip(&format!("{prefix}_2.fastq.gz")));
        let (mut q37, mut total, mut low, mut low_pairs) = (0u64, 0u64, 0u64, 0u64);
        for r in &reads {
            let q = r.qual.as_bytes();
            total += q.len() as u64;
            q37 += q.iter().filter(|&&c| c == b'F').count() as u64;
            for w in q.windows(2) {
                if w[0] == b',' {
                    low += 1;
                    low_pairs += (w[1] == b',') as u64;
                }
            }
        }
        (q37 as f64 / total as f64, low_pairs as f64 / low as f64)
    };
    let (q37_a, sticky_a) = summarise(&first);
    let (q37_b, sticky_b) = summarise(&second);
    assert!((q37_a - q37_b).abs() < 0.003, "{q37_a} vs {q37_b}");
    assert!(
        (sticky_a - sticky_b).abs() < 0.03,
        "{sticky_a} vs {sticky_b}"
    );
    assert!(sticky_b > 0.6);

    // --independent drops the dependence between neighbouring cycles
    let independent = scratch.path("independent.json");
    assert_success(&run(&[
        "profile",
        "-1",
        &format!("{first}_2.fastq.gz"),
        "-o",
        &independent,
        "--independent",
        "--max-reads",
        "10k",
        "--quiet",
    ]));
    let json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&independent).unwrap()).unwrap();
    assert!(json["reads"][0].get("transitions").is_none());
    assert_eq!(json["reads"].as_array().unwrap().len(), 1);
    let third = scratch.path("third");
    assert_success(&simulate(
        &third,
        &[
            "-n",
            "30000",
            "--read-length",
            "60",
            "--quality-profile",
            &independent,
        ],
    ));
    let (q37_c, sticky_c) = summarise(&third);
    assert!((q37_a - q37_c).abs() < 0.004, "{q37_a} vs {q37_c}");
    assert!(sticky_c < 0.1, "{sticky_c}");
}

#[test]
fn other_table_layouts_and_gzip_input() {
    let scratch = Scratch::new("tables");
    // a Salmon-style quant.sf with full GENCODE-style names, gzip-compressed
    let quant = scratch.path("quant.sf.gz");
    {
        use std::io::Write;
        let mut enc = flate2::write::GzEncoder::new(
            std::fs::File::create(&quant).unwrap(),
            flate2::Compression::default(),
        );
        writeln!(enc, "Name\tLength\tEffectiveLength\tTPM\tNumReads").unwrap();
        writeln!(
            enc,
            "TOY0003.1|TOYG002.1|-|-|BETA-201|BETA|900|protein_coding|\t900\t720\t10\t300"
        )
        .unwrap();
        writeln!(enc, "TOY0005.1\t350\t170\t30\t100").unwrap();
        enc.finish().unwrap();
    }
    let transcripts = data("transcripts.fa");
    let prefix = scratch.path("salmon");
    let out = run(&[
        "simulate",
        "-t",
        transcripts.to_str().unwrap(),
        "-a",
        &quant,
        "-o",
        &prefix,
        "-n",
        "1000",
        "--quiet",
    ]);
    assert_success(&out);
    let truth = std::fs::read_to_string(format!("{prefix}.truth.tsv")).unwrap();
    let counts: HashMap<&str, u64> = truth
        .lines()
        .skip(1)
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[0], f[5].parse().unwrap())
        })
        .collect();
    assert_eq!(counts["TOY0003.1"] + counts["TOY0005.1"], 1000);
    assert!(counts["TOY0003.1"] > counts["TOY0005.1"]);

    // the same table used as counts: 300 to 100
    let prefix = scratch.path("counts");
    let out = run(&[
        "simulate",
        "-t",
        transcripts.to_str().unwrap(),
        "-a",
        &quant,
        "-o",
        &prefix,
        "-n",
        "4000",
        "--abundance-column",
        "NumReads",
        "--abundance-type",
        "counts",
        "--quiet",
    ]);
    assert_success(&out);
    let truth = std::fs::read_to_string(format!("{prefix}.truth.tsv")).unwrap();
    let beta: u64 = truth
        .lines()
        .find(|l| l.starts_with("TOY0003.1"))
        .unwrap()
        .split('\t')
        .nth(5)
        .unwrap()
        .parse()
        .unwrap();
    assert!(
        (2850..3150).contains(&beta),
        "{beta} of 4000, expected about 3000"
    );
}

#[test]
fn errors_are_reported_without_a_panic() {
    let scratch = Scratch::new("errors");
    let prefix = scratch.path("x");
    let transcripts = data("transcripts.fa");
    let transcripts = transcripts.to_str().unwrap();

    let failure = |args: &[&str], code: i32, needle: &str| {
        let out = run(args);
        let stderr = String::from_utf8_lossy(&out.stderr).into_owned();
        assert_eq!(out.status.code(), Some(code), "args {args:?}: {stderr}");
        assert!(
            stderr.contains(needle),
            "args {args:?}: expected '{needle}' in: {stderr}"
        );
        assert!(!stderr.contains("panicked"), "{stderr}");
    };

    // identifiers that are not in the FASTA file
    let bad_table = scratch.path("bad.tsv");
    std::fs::write(
        &bad_table,
        "transcript_id\ttpm\nTOY0001.1\t1\nNOPE1\t2\nNOPE2\t3\n",
    )
    .unwrap();
    failure(
        &[
            "simulate",
            "-t",
            transcripts,
            "-a",
            &bad_table,
            "-o",
            &prefix,
        ],
        1,
        "2 transcript identifier(s) are not in the FASTA file: NOPE1, NOPE2",
    );
    // missing files
    let abundance = data("abundance.tsv");
    let abundance = abundance.to_str().unwrap();
    failure(
        &[
            "simulate",
            "-t",
            &scratch.path("missing.fa"),
            "-a",
            abundance,
            "-o",
            &prefix,
        ],
        1,
        "cannot open",
    );
    failure(
        &[
            "simulate",
            "-t",
            transcripts,
            "-a",
            &scratch.path("missing.tsv"),
            "-o",
            &prefix,
        ],
        1,
        "cannot open",
    );
    // output directory does not exist
    failure(
        &[
            "simulate",
            "-t",
            transcripts,
            "-a",
            abundance,
            "-o",
            &scratch.path("no/such/dir/x"),
        ],
        1,
        "cannot create",
    );
    // wrong column name
    failure(
        &[
            "simulate",
            "-t",
            transcripts,
            "-a",
            abundance,
            "-o",
            &prefix,
            "--abundance-column",
            "fpkm",
        ],
        1,
        "no abundance column 'fpkm'",
    );
    // invalid values
    failure(
        &[
            "simulate",
            "-t",
            transcripts,
            "-a",
            abundance,
            "-o",
            &prefix,
            "--fragment-min",
            "0",
        ],
        1,
        "minimum fragment length",
    );
    failure(
        &[
            "simulate",
            "-t",
            transcripts,
            "-a",
            abundance,
            "-o",
            &prefix,
            "--index1",
            "ACGU",
        ],
        1,
        "only A, C, G, T and N",
    );
    failure(
        &[
            "simulate",
            "-t",
            transcripts,
            "-a",
            abundance,
            "-o",
            &prefix,
            "--threads",
            "0",
        ],
        1,
        "at least 1",
    );
    failure(
        &[
            "simulate",
            "-t",
            transcripts,
            "-a",
            abundance,
            "-o",
            &prefix,
            "-q",
            abundance,
        ],
        1,
        "not a valid quality profile",
    );
    // a table in which nothing is expressed
    let zero_table = scratch.path("zero.tsv");
    std::fs::write(&zero_table, "transcript_id\ttpm\nTOY0001.1\t0\n").unwrap();
    failure(
        &[
            "simulate",
            "-t",
            transcripts,
            "-a",
            &zero_table,
            "-o",
            &prefix,
        ],
        1,
        "nothing to simulate",
    );
    // command-line mistakes: status 2
    failure(&["simulate", "-t", transcripts], 2, "required");
    failure(
        &[
            "simulate",
            "-t",
            transcripts,
            "-a",
            abundance,
            "-o",
            &prefix,
            "--library-type",
            "sideways",
        ],
        2,
        "invalid value",
    );
    failure(
        &[
            "simulate",
            "-t",
            transcripts,
            "-a",
            abundance,
            "-o",
            &prefix,
            "--compression-level",
            "12",
        ],
        2,
        "invalid value",
    );
    failure(&["frobnicate"], 2, "unrecognized subcommand");
    // profile: not a FASTQ file
    failure(
        &["profile", "-1", transcripts, "-o", &scratch.path("p.json")],
        1,
        "expected a FASTQ header",
    );

    // --help and --version work
    let help = run(&["simulate", "--help"]);
    assert!(help.status.success());
    let text = String::from_utf8_lossy(&help.stdout);
    for option in [
        "--library-type",
        "--fragment-mean",
        "--adapter-fill",
        "--quality-profile",
        "--seed",
    ] {
        assert!(text.contains(option), "help does not mention {option}");
    }
    let version = run(&["--version"]);
    assert!(String::from_utf8_lossy(&version.stdout).contains(env!("CARGO_PKG_VERSION")));
}
