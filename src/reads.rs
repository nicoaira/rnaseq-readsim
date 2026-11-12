//! Turning a fragment into the two error-free reads: strand, reverse
//! complement and adapter read-through.

use crate::rng::Rng;

/// Default sequence that read 1 runs into past the 3′ end of a short
/// fragment: the Illumina TruSeq index (i7) adapter, with `{i7}` standing for
/// the index. Source: Illumina Adapter Sequences, document 1000000002694
/// (TruSeq single and unique-dual indexes).
pub const DEFAULT_ADAPTER1: &str = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC{i7}ATCTCGTATGCCGTCTTCTGCTTG";

/// Default sequence that read 2 runs into: the reverse complement of the
/// Illumina TruSeq unique-dual-index i5 adapter, with `{i5}` standing for the
/// index as read 2 reads it. Same source as [`DEFAULT_ADAPTER1`].
pub const DEFAULT_ADAPTER2: &str =
    "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT{i5}GTGTAGATCTCGGTGGTCGCCGTATCATT";

/// What read 2 runs into in a single-index TruSeq library: the reverse
/// complement of the TruSeq universal adapter (no index). Not the default;
/// pass it as `--adapter2` together with an empty `--index2`.
pub const TRUSEQ_UNIVERSAL_READ2: &str =
    "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT";

/// Default i7 index (IDT for Illumina TruSeq UD index UDI0001).
pub const DEFAULT_INDEX1: &str = "CCGCGGTT";

/// Default i5 index (UDI0001, in the orientation in which read 2 reads it).
pub const DEFAULT_INDEX2: &str = "CTAGCGCT";

/// Which strand the reads come from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LibraryType {
    /// dUTP / "fr-firststrand" libraries: read 1 is the reverse complement of
    /// the transcript, read 2 is on the transcript strand (Salmon: ISR).
    Reverse,
    /// Read 1 is on the transcript strand, read 2 is the reverse complement
    /// (Salmon: ISF).
    Forward,
    /// Either orientation with probability one half (Salmon: IU).
    Unstranded,
}

impl LibraryType {
    pub fn name(self) -> &'static str {
        match self {
            LibraryType::Reverse => "reverse",
            LibraryType::Forward => "forward",
            LibraryType::Unstranded => "unstranded",
        }
    }

    /// Whether read 1 of the next fragment is on the transcript strand.
    /// Only unstranded libraries use a random number.
    #[inline]
    pub fn read1_on_transcript_strand(self, rng: &mut Rng) -> bool {
        match self {
            LibraryType::Forward => true,
            LibraryType::Reverse => false,
            LibraryType::Unstranded => rng.next_u64() >> 63 == 1,
        }
    }
}

/// What a read contains once it has run past the end of the adapter.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AdapterFill {
    /// `G`, the "no signal" call of two-colour instruments.
    PolyG,
    /// `A`.
    PolyA,
    /// Uniformly random bases.
    Random,
}

impl AdapterFill {
    pub fn name(self) -> &'static str {
        match self {
            AdapterFill::PolyG => "polyG",
            AdapterFill::PolyA => "polyA",
            AdapterFill::Random => "random",
        }
    }
}

const fn complement_table() -> [u8; 256] {
    let mut t = [b'N'; 256];
    t[b'A' as usize] = b'T';
    t[b'C' as usize] = b'G';
    t[b'G' as usize] = b'C';
    t[b'T' as usize] = b'A';
    t
}

/// Complement of an upper-case base; anything that is not `ACGT` gives `N`.
pub static COMPLEMENT: [u8; 256] = complement_table();

/// The reverse complement of an upper-case sequence.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| COMPLEMENT[b as usize]).collect()
}

/// Checks and normalises a nucleotide sequence given on the command line
/// (adapter or index): upper case, only `ACGTN`.
pub fn parse_nucleotides(text: &str, what: &str) -> Result<Vec<u8>, String> {
    let mut out = Vec::with_capacity(text.len());
    for c in text.bytes() {
        let u = c.to_ascii_uppercase();
        if !matches!(u, b'A' | b'C' | b'G' | b'T' | b'N') {
            return Err(format!(
                "{what} '{text}' contains '{}'; only A, C, G, T and N are allowed",
                (c as char).escape_default()
            ));
        }
        out.push(u);
    }
    Ok(out)
}

/// Replaces the placeholders `{i7}` and `{i5}` of an adapter template by the
/// index sequences and checks the result.
pub fn resolve_adapter(
    template: &str,
    index1: &[u8],
    index2: &[u8],
    what: &str,
) -> Result<Vec<u8>, String> {
    let i7 = String::from_utf8_lossy(index1);
    let i5 = String::from_utf8_lossy(index2);
    let filled = template.replace("{i7}", &i7).replace("{i5}", &i5);
    if filled.contains('{') || filled.contains('}') {
        return Err(format!(
            "{what} '{template}' contains a placeholder other than {{i7}} or {{i5}}"
        ));
    }
    parse_nucleotides(&filled, what)
}

/// Writes one error-free read into `out` (whose length is the read length).
///
/// `fragment` is given on the transcript strand. A read on the transcript
/// strand (`on_transcript_strand`) starts at the 5′ end of the fragment; the
/// other read is the reverse complement and starts at its 3′ end. If the
/// fragment is shorter than the read, the read continues into `adapter` and,
/// past its end, into the filler.
#[inline]
pub fn build_read(
    out: &mut [u8],
    fragment: &[u8],
    on_transcript_strand: bool,
    adapter: &[u8],
    fill: AdapterFill,
    rng: &mut Rng,
) {
    let read_len = out.len();
    let n = fragment.len().min(read_len);
    if on_transcript_strand {
        out[..n].copy_from_slice(&fragment[..n]);
    } else {
        for (o, &b) in out[..n].iter_mut().zip(fragment.iter().rev()) {
            *o = COMPLEMENT[b as usize];
        }
    }
    if n < read_len {
        let a = adapter.len().min(read_len - n);
        out[n..n + a].copy_from_slice(&adapter[..a]);
        let rest = &mut out[n + a..];
        match fill {
            AdapterFill::PolyG => rest.fill(b'G'),
            AdapterFill::PolyA => rest.fill(b'A'),
            AdapterFill::Random => {
                for o in rest {
                    *o = b"ACGT"[(rng.next_u64() >> 62) as usize];
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reverse_complement_basics() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AACCGGTT"), b"AACCGGTT");
        assert_eq!(reverse_complement(b"AAAC"), b"GTTT");
        assert_eq!(reverse_complement(b"ANGT"), b"ACNT");
        assert_eq!(reverse_complement(b""), b"");
        let s = b"GATTACAGATTACA";
        assert_eq!(reverse_complement(&reverse_complement(s)), s);
    }

    #[test]
    fn default_adapters_match_the_published_sequences() {
        // Illumina Adapter Sequences (document 1000000002694): the trimming
        // sequences for TruSeq libraries are prefixes of the read-through.
        assert!(DEFAULT_ADAPTER1.starts_with("AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"));
        assert!(DEFAULT_ADAPTER2.starts_with("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"));
        // Read 1 runs into the index adapter
        // (A)GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[i7]ATCTCGTATGCCGTCTTCTGCTTG.
        let a1 = resolve_adapter(DEFAULT_ADAPTER1, b"ATCACG", b"", "adapter 1").unwrap();
        assert_eq!(
            a1,
            b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"
        );
        // Read 2 of a unique-dual-index library runs into the reverse
        // complement of the i5 adapter
        // AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT.
        // With the unique dual index UDI0001 (i5 in the adapter: AGCGCTAG).
        let dual = b"AATGATACGGCGACCACCGAGATCTACACAGCGCTAGACACTCTTTCCCTACACGACGCTCTTCCGATCT";
        let a2 = resolve_adapter(
            DEFAULT_ADAPTER2,
            DEFAULT_INDEX1.as_bytes(),
            DEFAULT_INDEX2.as_bytes(),
            "adapter 2",
        )
        .unwrap();
        assert_eq!(a2, reverse_complement(dual));
        // Single-index TruSeq libraries carry the universal adapter, which
        // has no index and no "ACAC" spacer; its reverse complement can be
        // given as it is.
        let universal = b"AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
        let a2 = resolve_adapter(TRUSEQ_UNIVERSAL_READ2, b"ATCACG", b"", "adapter 2").unwrap();
        assert_eq!(a2, reverse_complement(universal));
        assert!(TRUSEQ_UNIVERSAL_READ2.starts_with("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"));
    }

    #[test]
    fn adapter_templates() {
        assert_eq!(
            resolve_adapter("ac{i7}gt{i5}", b"TT", b"GG", "adapter").unwrap(),
            b"ACTTGTGG"
        );
        assert_eq!(resolve_adapter("", b"TT", b"GG", "adapter").unwrap(), b"");
        assert!(resolve_adapter("AC{i9}", b"TT", b"GG", "adapter").is_err());
        assert!(resolve_adapter("ACXT", b"TT", b"GG", "adapter").is_err());
        assert!(parse_nucleotides("ACGU", "index").is_err());
        assert_eq!(parse_nucleotides("acgtn", "index").unwrap(), b"ACGTN");
    }

    #[test]
    fn strand_per_library_type() {
        let mut rng = Rng::for_stream(1, 0);
        for _ in 0..100 {
            assert!(LibraryType::Forward.read1_on_transcript_strand(&mut rng));
            assert!(!LibraryType::Reverse.read1_on_transcript_strand(&mut rng));
        }
        let n = 100_000;
        let forward = (0..n)
            .filter(|_| LibraryType::Unstranded.read1_on_transcript_strand(&mut rng))
            .count();
        assert!((49_000..51_000).contains(&forward), "{forward}");
    }

    #[test]
    fn reads_of_a_long_fragment() {
        let fragment = b"AAAACCCCGGGGTTTTACGT"; // 20 bases
        let mut rng = Rng::for_stream(1, 0);
        let mut read = [0u8; 8];
        build_read(
            &mut read,
            fragment,
            true,
            b"TTTT",
            AdapterFill::PolyG,
            &mut rng,
        );
        assert_eq!(&read, b"AAAACCCC");
        build_read(
            &mut read,
            fragment,
            false,
            b"TTTT",
            AdapterFill::PolyG,
            &mut rng,
        );
        // reverse complement of the last 8 bases TTTTACGT
        assert_eq!(&read, b"ACGTAAAA");
    }

    #[test]
    fn read_through_geometry() {
        let fragment = b"AACCG"; // 5 bases, read length 12
        let adapter = b"TTTT";
        let mut rng = Rng::for_stream(1, 0);
        let mut read = [0u8; 12];
        build_read(
            &mut read,
            fragment,
            true,
            adapter,
            AdapterFill::PolyG,
            &mut rng,
        );
        assert_eq!(&read, b"AACCGTTTTGGG");
        build_read(
            &mut read,
            fragment,
            false,
            adapter,
            AdapterFill::PolyA,
            &mut rng,
        );
        assert_eq!(&read, b"CGGTTTTTTAAA");
        // fragment exactly as long as the read: no adapter
        let mut read = [0u8; 5];
        build_read(
            &mut read,
            fragment,
            true,
            adapter,
            AdapterFill::PolyG,
            &mut rng,
        );
        assert_eq!(&read, b"AACCG");
        // one base shorter than the read: one adapter base
        let mut read = [0u8; 6];
        build_read(
            &mut read,
            fragment,
            true,
            adapter,
            AdapterFill::PolyG,
            &mut rng,
        );
        assert_eq!(&read, b"AACCGT");
        // adapter longer than the rest of the read: it is cut
        let mut read = [0u8; 7];
        build_read(
            &mut read,
            fragment,
            false,
            b"ACGTACGT",
            AdapterFill::PolyG,
            &mut rng,
        );
        assert_eq!(&read, b"CGGTTAC");
        // empty adapter: only filler
        let mut read = [0u8; 8];
        build_read(&mut read, fragment, true, b"", AdapterFill::PolyA, &mut rng);
        assert_eq!(&read, b"AACCGAAA");
        // random filler: valid bases, insert and adapter untouched
        let mut read = [0u8; 40];
        build_read(
            &mut read,
            fragment,
            true,
            adapter,
            AdapterFill::Random,
            &mut rng,
        );
        assert_eq!(&read[..9], b"AACCGTTTT");
        assert!(read[9..].iter().all(|b| b"ACGT".contains(b)));
        let distinct: std::collections::HashSet<u8> = read[9..].iter().copied().collect();
        assert!(distinct.len() > 1);
    }

    #[test]
    fn mates_of_a_short_fragment_are_reverse_complements_over_the_insert() {
        let fragment = b"ACGGTCA";
        let mut rng = Rng::for_stream(1, 0);
        let mut r1 = [0u8; 10];
        let mut r2 = [0u8; 10];
        build_read(
            &mut r1,
            fragment,
            false,
            b"AGATCGG",
            AdapterFill::PolyG,
            &mut rng,
        );
        build_read(
            &mut r2,
            fragment,
            true,
            b"AGATCGG",
            AdapterFill::PolyG,
            &mut rng,
        );
        assert_eq!(&r2[..7], fragment);
        assert_eq!(&r1[..7], reverse_complement(fragment).as_slice());
        assert_eq!(&r1[7..], b"AGA");
        assert_eq!(&r2[7..], b"AGA");
    }
}
