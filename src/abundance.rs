//! Reading the abundance table (tab-separated, with a header line).

use std::io::BufRead;
use std::path::Path;

use crate::error::{Context, Error, Result};
use crate::fasta::{Transcriptome, open_text, parse_id};

/// Identifier columns tried, in this order, when none is named.
pub const DEFAULT_ID_COLUMNS: [&str; 3] = ["transcript_id", "Name", "target_id"];
/// Abundance columns tried when none is named (matched without regard to case).
pub const DEFAULT_ABUNDANCE_COLUMNS: [&str; 1] = ["tpm"];

/// The abundances of all transcripts of a transcriptome.
#[derive(Debug, Clone)]
pub struct Abundances {
    /// One value per transcript, in FASTA order; 0 for transcripts that are
    /// not in the table.
    pub values: Vec<f64>,
    /// Number of data rows in the table.
    pub rows: usize,
    /// Name of the identifier column that was used.
    pub id_column: String,
    /// Name of the abundance column that was used.
    pub abundance_column: String,
}

fn find_column(header: &[&str], wanted: &[&str], what: &str, explicit: bool) -> Result<usize> {
    for w in wanted {
        if let Some(i) = header.iter().position(|h| h == w) {
            return Ok(i);
        }
    }
    for w in wanted {
        let hits: Vec<usize> = header
            .iter()
            .enumerate()
            .filter(|(_, h)| h.eq_ignore_ascii_case(w))
            .map(|(i, _)| i)
            .collect();
        match hits.len() {
            0 => {}
            1 => return Ok(hits[0]),
            _ => {
                return Err(Error(format!(
                    "the abundance table has several columns named '{w}' when case is ignored; \
                     name the {what} column exactly"
                )));
            }
        }
    }
    let hint = if header.len() == 1 {
        " (only one column was found: is the table tab-separated?)"
    } else {
        ""
    };
    let wanted_text = wanted
        .iter()
        .map(|w| format!("'{w}'"))
        .collect::<Vec<_>>()
        .join(" or ");
    let option = if explicit {
        String::new()
    } else {
        format!(
            "; name the column with --{}-column",
            if what == "identifier" {
                "id"
            } else {
                "abundance"
            }
        )
    };
    Err(Error(format!(
        "the abundance table has no {what} column {wanted_text}{hint}; its columns are: {}{option}",
        header.join(", ")
    )))
}

/// Parses an abundance table from a reader.
///
/// `name` is used in error messages. Empty lines and lines starting with `#`
/// are skipped. Identifiers are cut at the first white-space character or
/// `|`, like FASTA headers.
pub fn read_abundances_from<R: BufRead>(
    mut reader: R,
    name: &str,
    transcriptome: &Transcriptome,
    id_column: Option<&str>,
    abundance_column: Option<&str>,
) -> Result<Abundances> {
    let index = transcriptome.index()?;
    let mut values = vec![0.0f64; transcriptome.len()];
    let mut seen = vec![false; transcriptome.len()];
    let mut columns: Option<(usize, usize, String, String)> = None;
    let mut unknown: Vec<String> = Vec::new();
    let mut unknown_count = 0usize;
    let mut rows = 0usize;
    let mut line = String::new();
    let mut line_no = 0u64;
    loop {
        line.clear();
        let n = reader
            .read_line(&mut line)
            .context(|| format!("cannot read {name}"))?;
        if n == 0 {
            break;
        }
        line_no += 1;
        let text = line.trim_end_matches(['\n', '\r']);
        if text.trim().is_empty() || text.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = text.split('\t').collect();
        let Some((id_col, ab_col, _, _)) = &columns else {
            let header: Vec<&str> = fields.iter().map(|f| f.trim()).collect();
            let id_col = match id_column {
                Some(c) => find_column(&header, &[c], "identifier", true)?,
                None => find_column(&header, &DEFAULT_ID_COLUMNS, "identifier", false)?,
            };
            let ab_col = match abundance_column {
                Some(c) => find_column(&header, &[c], "abundance", true)?,
                None => find_column(&header, &DEFAULT_ABUNDANCE_COLUMNS, "abundance", false)?,
            };
            if id_col == ab_col {
                return Err(Error(format!(
                    "{name}: the identifier column and the abundance column are the same column"
                )));
            }
            columns = Some((
                id_col,
                ab_col,
                header[id_col].to_owned(),
                header[ab_col].to_owned(),
            ));
            continue;
        };
        let needed = (*id_col).max(*ab_col);
        if fields.len() <= needed {
            return Err(Error(format!(
                "{name}: line {line_no} has {} columns, fewer than the header",
                fields.len()
            )));
        }
        rows += 1;
        let raw_id = fields[*id_col].trim();
        let id = std::str::from_utf8(parse_id(raw_id.as_bytes())).unwrap_or(raw_id);
        if id.is_empty() {
            return Err(Error(format!(
                "{name}: line {line_no}: empty transcript identifier"
            )));
        }
        let raw_value = fields[*ab_col].trim();
        let value: f64 = raw_value.parse().map_err(|_| {
            Error(format!(
                "{name}: line {line_no}: abundance '{raw_value}' of '{id}' is not a number"
            ))
        })?;
        if !value.is_finite() || value < 0.0 {
            return Err(Error(format!(
                "{name}: line {line_no}: abundance '{raw_value}' of '{id}' is not a finite, \
                 non-negative number"
            )));
        }
        match index.get(id) {
            Some(&i) => {
                if seen[i] {
                    return Err(Error(format!(
                        "{name}: line {line_no}: transcript '{id}' occurs more than once"
                    )));
                }
                seen[i] = true;
                values[i] = value;
            }
            None => {
                unknown_count += 1;
                if unknown.len() < 5 {
                    unknown.push(id.to_owned());
                }
            }
        }
    }
    let Some((_, _, id_name, abundance_name)) = columns else {
        return Err(Error(format!("{name}: the table is empty")));
    };
    if unknown_count > 0 {
        let more = if unknown_count > unknown.len() {
            ", ..."
        } else {
            ""
        };
        return Err(Error(format!(
            "{name}: {unknown_count} transcript identifier(s) are not in the FASTA file: {}{more}",
            unknown.join(", ")
        )));
    }
    Ok(Abundances {
        values,
        rows,
        id_column: id_name,
        abundance_column: abundance_name,
    })
}

/// Reads an abundance table from a file (plain or gzip).
pub fn read_abundances(
    path: &Path,
    transcriptome: &Transcriptome,
    id_column: Option<&str>,
    abundance_column: Option<&str>,
) -> Result<Abundances> {
    let reader = open_text(path)?;
    read_abundances_from(
        reader,
        &path.display().to_string(),
        transcriptome,
        id_column,
        abundance_column,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tx() -> Transcriptome {
        Transcriptome::from_records([("t1", "ACGT"), ("t2", "ACGT"), ("t3", "ACGT")]).unwrap()
    }

    fn parse(text: &str, id: Option<&str>, ab: Option<&str>) -> Result<Abundances> {
        read_abundances_from(text.as_bytes(), "table", &tx(), id, ab)
    }

    #[test]
    fn default_columns() {
        let a = parse("transcript_id\ttpm\nt1\t5\nt3\t2.5\n", None, None).unwrap();
        assert_eq!(a.values, vec![5.0, 0.0, 2.5]);
        assert_eq!(a.rows, 2);
        assert_eq!(a.id_column, "transcript_id");
    }

    #[test]
    fn salmon_rsem_and_kallisto_layouts() {
        let salmon = "Name\tLength\tEffectiveLength\tTPM\tNumReads\nt2\t4\t4\t7.5\t3\n";
        let a = parse(salmon, None, None).unwrap();
        assert_eq!(a.values, vec![0.0, 7.5, 0.0]);
        assert_eq!(a.abundance_column, "TPM");
        let rsem = "transcript_id\tgene_id\tlength\teffective_length\texpected_count\tTPM\tFPKM\tIsoPct\nt1\tg\t4\t4\t9\t1e2\t1\t100\n";
        let a = parse(rsem, None, None).unwrap();
        assert_eq!(a.values, vec![100.0, 0.0, 0.0]);
        let a = parse(rsem, None, Some("expected_count")).unwrap();
        assert_eq!(a.values, vec![9.0, 0.0, 0.0]);
        let kallisto = "target_id\tlength\teff_length\test_counts\ttpm\nt3\t4\t4\t1\t0.25\n";
        let a = parse(kallisto, None, None).unwrap();
        assert_eq!(a.values, vec![0.0, 0.0, 0.25]);
    }

    #[test]
    fn identifiers_are_cut_like_fasta_headers() {
        let a = parse("Name\tTPM\nt1|gene|x\t3\n", None, None).unwrap();
        assert_eq!(a.values[0], 3.0);
    }

    #[test]
    fn comments_blank_lines_and_crlf() {
        let a = parse(
            "# made by hand\r\n\r\ntranscript_id\ttpm\r\nt2\t1\r\n",
            None,
            None,
        )
        .unwrap();
        assert_eq!(a.values, vec![0.0, 1.0, 0.0]);
    }

    #[test]
    fn errors() {
        let unknown = parse("transcript_id\ttpm\nx1\t1\nx2\t1\nt1\t1\n", None, None).unwrap_err();
        assert!(
            unknown.0.contains("2 transcript identifier(s)"),
            "{unknown}"
        );
        assert!(unknown.0.contains("x1, x2"), "{unknown}");
        assert!(parse("transcript_id\ttpm\nt1\t1\nt1\t2\n", None, None).is_err());
        assert!(parse("transcript_id\ttpm\nt1\tabc\n", None, None).is_err());
        assert!(parse("transcript_id\ttpm\nt1\t-1\n", None, None).is_err());
        assert!(parse("transcript_id\ttpm\nt1\tNaN\n", None, None).is_err());
        assert!(parse("transcript_id\ttpm\nt1\n", None, None).is_err());
        assert!(parse("", None, None).is_err());
        let missing = parse("id,tpm\nt1,1\n", None, None).unwrap_err();
        assert!(missing.0.contains("tab-separated"), "{missing}");
        let missing = parse("a\tb\nt1\t1\n", Some("zzz"), None).unwrap_err();
        assert!(missing.0.contains("'zzz'"), "{missing}");
    }

    #[test]
    fn many_unknown_identifiers_are_summarised() {
        let mut text = String::from("transcript_id\ttpm\n");
        for i in 0..9 {
            text.push_str(&format!("u{i}\t1\n"));
        }
        let e = parse(&text, None, None).unwrap_err();
        assert!(e.0.contains("9 transcript identifier(s)"), "{e}");
        assert!(e.0.contains("u0, u1, u2, u3, u4, ..."), "{e}");
    }
}
