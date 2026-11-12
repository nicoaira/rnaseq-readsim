//! Reading the transcript FASTA (plain or gzip) into memory.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::bufread::MultiGzDecoder;

use crate::error::{Context, Error, Result};

const READ_BUFFER: usize = 1 << 20;

/// Opens a file for buffered reading; gzip is recognised by its magic bytes,
/// not by the file name.
pub fn open_text(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).context(|| format!("cannot open {}", path.display()))?;
    let mut reader = BufReader::with_capacity(READ_BUFFER, file);
    let is_gzip = {
        let head = reader
            .fill_buf()
            .context(|| format!("cannot read {}", path.display()))?;
        head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b
    };
    if is_gzip {
        Ok(Box::new(BufReader::with_capacity(
            READ_BUFFER,
            MultiGzDecoder::new(reader),
        )))
    } else {
        Ok(Box::new(reader))
    }
}

/// The sequence identifier of a header or table field: everything up to the
/// first white-space character or `|`.
///
/// `ENST00000832824.1|ENSG00000290825.2|-|...` and `ENST00000832824.1 extra`
/// both give `ENST00000832824.1`.
pub fn parse_id(field: &[u8]) -> &[u8] {
    let end = field
        .iter()
        .position(|&b| b == b'|' || b.is_ascii_whitespace())
        .unwrap_or(field.len());
    &field[..end]
}

/// Maps a FASTA sequence character to what is stored: `ACGT` in upper case,
/// `U` as `T`, any other IUPAC nucleotide code as `N`, and 0 for a character
/// that is not a nucleotide code.
const fn normalisation_table() -> [u8; 256] {
    let mut t = [0u8; 256];
    let plain = b"ACGT";
    let mut i = 0;
    while i < plain.len() {
        t[plain[i] as usize] = plain[i];
        t[plain[i].to_ascii_lowercase() as usize] = plain[i];
        i += 1;
    }
    t[b'U' as usize] = b'T';
    t[b'u' as usize] = b'T';
    let ambiguous = b"NRYSWKMBDHV";
    let mut i = 0;
    while i < ambiguous.len() {
        t[ambiguous[i] as usize] = b'N';
        t[ambiguous[i].to_ascii_lowercase() as usize] = b'N';
        i += 1;
    }
    t
}

static NORMALISE: [u8; 256] = normalisation_table();

/// All transcripts, concatenated in one buffer.
#[derive(Debug, Default)]
pub struct Transcriptome {
    ids: Vec<String>,
    seq: Vec<u8>,
    /// `offsets[i]..offsets[i + 1]` is transcript `i` in `seq`.
    offsets: Vec<usize>,
}

impl Transcriptome {
    /// Number of transcripts.
    pub fn len(&self) -> usize {
        self.ids.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }

    pub fn id(&self, i: usize) -> &str {
        &self.ids[i]
    }

    pub fn ids(&self) -> &[String] {
        &self.ids
    }

    #[inline]
    pub fn sequence(&self, i: usize) -> &[u8] {
        &self.seq[self.offsets[i]..self.offsets[i + 1]]
    }

    #[inline]
    pub fn length(&self, i: usize) -> usize {
        self.offsets[i + 1] - self.offsets[i]
    }

    pub fn total_bases(&self) -> usize {
        self.seq.len()
    }

    /// Builds a transcriptome from `(id, sequence)` pairs (used by tests and
    /// by callers that already hold sequences in memory).
    pub fn from_records<I, S, B>(records: I) -> Result<Self>
    where
        I: IntoIterator<Item = (S, B)>,
        S: Into<String>,
        B: AsRef<[u8]>,
    {
        let mut t = Transcriptome {
            offsets: vec![0],
            ..Default::default()
        };
        for (id, seq) in records {
            let id: String = id.into();
            if seq.as_ref().is_empty() {
                return Err(Error(format!("transcript '{id}' has an empty sequence")));
            }
            for &b in seq.as_ref() {
                let n = NORMALISE[b as usize];
                if n == 0 {
                    return Err(Error(format!(
                        "transcript '{id}' contains the character '{}', which is not a nucleotide code",
                        (b as char).escape_default()
                    )));
                }
                t.seq.push(n);
            }
            t.ids.push(id);
            t.offsets.push(t.seq.len());
        }
        t.check()?;
        Ok(t)
    }

    /// Reads a FASTA file (plain or gzip).
    pub fn from_fasta(path: &Path) -> Result<Self> {
        let mut reader = open_text(path)?;
        let mut t = Transcriptome {
            offsets: vec![0],
            ..Default::default()
        };
        if let Ok(meta) = std::fs::metadata(path) {
            // A plain FASTA is slightly larger than its sequence; for gzip
            // this is only a first guess.
            t.seq.reserve(meta.len() as usize);
        }
        let mut line: Vec<u8> = Vec::with_capacity(1 << 16);
        let mut line_no = 0u64;
        let mut in_record = false;
        loop {
            line.clear();
            let n = reader
                .read_until(b'\n', &mut line)
                .context(|| format!("cannot read {}", path.display()))?;
            if n == 0 {
                break;
            }
            line_no += 1;
            while matches!(line.last(), Some(b'\n' | b'\r' | b' ' | b'\t')) {
                line.pop();
            }
            if line.is_empty() {
                continue;
            }
            if line[0] == b'>' {
                if in_record {
                    t.finish_record(path)?;
                }
                let id = parse_id(&line[1..]);
                if id.is_empty() {
                    return Err(Error(format!(
                        "{}: line {line_no}: the header has no sequence identifier",
                        path.display()
                    )));
                }
                let id = std::str::from_utf8(id).map_err(|_| {
                    Error(format!(
                        "{}: line {line_no}: the sequence identifier is not valid UTF-8",
                        path.display()
                    ))
                })?;
                t.ids.push(id.to_owned());
                in_record = true;
            } else if line[0] == b';' {
                // old-style FASTA comment line
                continue;
            } else {
                if !in_record {
                    return Err(Error(format!(
                        "{}: line {line_no}: expected a FASTA header starting with '>'",
                        path.display()
                    )));
                }
                let start = t.seq.len();
                t.seq.extend(line.iter().map(|&b| NORMALISE[b as usize]));
                if t.seq[start..].contains(&0) {
                    let bad = line
                        .iter()
                        .find(|&&b| NORMALISE[b as usize] == 0)
                        .copied()
                        .unwrap_or(b'?');
                    return Err(Error(format!(
                        "{}: line {line_no}: '{}' is not a nucleotide code (sequence '{}')",
                        path.display(),
                        (bad as char).escape_default(),
                        t.ids.last().map(String::as_str).unwrap_or("")
                    )));
                }
            }
        }
        if in_record {
            t.finish_record(path)?;
        }
        if t.ids.is_empty() {
            return Err(Error(format!("{}: no FASTA records found", path.display())));
        }
        t.seq.shrink_to_fit();
        t.check()?;
        Ok(t)
    }

    fn finish_record(&mut self, path: &Path) -> Result<()> {
        let start = *self.offsets.last().expect("offsets start with 0");
        if self.seq.len() == start {
            return Err(Error(format!(
                "{}: sequence '{}' is empty",
                path.display(),
                self.ids.last().map(String::as_str).unwrap_or("")
            )));
        }
        self.offsets.push(self.seq.len());
        Ok(())
    }

    fn check(&self) -> Result<()> {
        if self.ids.len() > u32::MAX as usize {
            return Err(Error("too many sequences (more than 2^32 - 1)".into()));
        }
        self.index().map(|_| ())
    }

    /// Maps identifiers to transcript numbers; fails on a duplicated
    /// identifier.
    pub fn index(&self) -> Result<HashMap<&str, usize>> {
        let mut map = HashMap::with_capacity(self.ids.len());
        for (i, id) in self.ids.iter().enumerate() {
            if map.insert(id.as_str(), i).is_some() {
                return Err(Error(format!(
                    "the sequence identifier '{id}' occurs more than once in the FASTA file"
                )));
            }
        }
        Ok(map)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn id_parsing() {
        assert_eq!(
            parse_id(b"ENST00000832824.1|ENSG00000290825.2|-|-|DDX11L16-260|DDX11L16|1379|lncRNA|"),
            b"ENST00000832824.1"
        );
        assert_eq!(parse_id(b"ENST00000832824.1"), b"ENST00000832824.1");
        assert_eq!(parse_id(b"tx1 gene=abc"), b"tx1");
        assert_eq!(parse_id(b"tx1\tgene=abc"), b"tx1");
        assert_eq!(parse_id(b"tx1 a|b"), b"tx1");
        assert_eq!(parse_id(b""), b"");
        assert_eq!(parse_id(b"|x"), b"");
    }

    fn temp_path(name: &str) -> std::path::PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "rnaseq-readsim-fasta-{}-{name}",
            std::process::id()
        ));
        std::fs::create_dir_all(&dir).unwrap();
        dir.join("in.fa")
    }

    #[test]
    fn reads_plain_and_gzip() {
        let text = b">t1|gene1|x\nACGTac\ngtNN\r\n\n>t2 description\nuuRY\n";
        let plain = temp_path("plain");
        std::fs::write(&plain, text).unwrap();
        let gz = temp_path("gz");
        let mut enc = flate2::write::GzEncoder::new(
            File::create(&gz).unwrap(),
            flate2::Compression::default(),
        );
        enc.write_all(text).unwrap();
        enc.finish().unwrap();
        for path in [&plain, &gz] {
            let t = Transcriptome::from_fasta(path).unwrap();
            assert_eq!(t.len(), 2);
            assert_eq!(t.id(0), "t1");
            assert_eq!(t.id(1), "t2");
            assert_eq!(t.sequence(0), b"ACGTACGTNN");
            assert_eq!(t.sequence(1), b"TTNN");
            assert_eq!(t.length(0), 10);
            assert_eq!(t.total_bases(), 14);
        }
        std::fs::remove_dir_all(plain.parent().unwrap()).ok();
        std::fs::remove_dir_all(gz.parent().unwrap()).ok();
    }

    #[test]
    fn rejects_malformed_input() {
        let cases: [(&str, &[u8]); 5] = [
            ("nohdr", b"ACGT\n"),
            ("dup", b">a\nAC\n>a\nGT\n"),
            ("empty", b">a\n>b\nAC\n"),
            ("badchar", b">a\nAC1T\n"),
            ("nothing", b"\n\n"),
        ];
        for (name, text) in cases {
            let path = temp_path(name);
            std::fs::write(&path, text).unwrap();
            assert!(Transcriptome::from_fasta(&path).is_err(), "case {name}");
            std::fs::remove_dir_all(path.parent().unwrap()).ok();
        }
        assert!(Transcriptome::from_fasta(Path::new("/nonexistent/file.fa")).is_err());
    }

    #[test]
    fn from_records_checks_input() {
        let t = Transcriptome::from_records([("a", "acgu"), ("b", "NNN")]).unwrap();
        assert_eq!(t.sequence(0), b"ACGT");
        assert!(Transcriptome::from_records([("a", "AC"), ("a", "GT")]).is_err());
        assert!(Transcriptome::from_records([("a", "")]).is_err());
        assert!(Transcriptome::from_records([("a", "AC-T")]).is_err());
    }
}
