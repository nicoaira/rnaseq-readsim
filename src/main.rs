//! Command-line interface of rnaseq-readsim.

use std::path::PathBuf;
use std::process::ExitCode;

use clap::{Args, Parser, Subcommand, ValueEnum};

use rnaseq_readsim::error::{Context, Error, Result};
use rnaseq_readsim::quality::{
    PROFILE_FORMAT, PROFILE_VERSION, Profile, ReadProfile, learn_from_fastq,
};
use rnaseq_readsim::reads::{
    AdapterFill, DEFAULT_ADAPTER1, DEFAULT_ADAPTER2, DEFAULT_INDEX1, DEFAULT_INDEX2, LibraryType,
};
use rnaseq_readsim::simulate::{self, AbundanceType, Config, ReadParameters};

/// Simulate Illumina paired-end bulk RNA-seq reads from a transcriptome and
/// an abundance table.
///
/// Exit status: 0 on success, 1 if the run failed (unreadable or invalid
/// input, output that cannot be written), 2 if the command line is invalid.
#[derive(Parser, Debug)]
#[command(
    name = "rnaseq-readsim",
    version,
    about,
    long_about,
    max_term_width = 100
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Simulate paired-end reads
    ///
    /// Writes <PREFIX>_1.fastq.gz and <PREFIX>_2.fastq.gz, the per-transcript
    /// truth table <PREFIX>.truth.tsv and the run report <PREFIX>.report.json.
    /// The same inputs, options and seed give byte-identical files, whatever
    /// the number of threads.
    Simulate(Box<SimulateArgs>),
    /// Learn a quality profile from real FASTQ files
    ///
    /// Counts, for every cycle, how often each quality value occurs and
    /// (unless --independent is given) how often each value follows each
    /// value of the previous cycle. The result is a small JSON file for
    /// `simulate --quality-profile`.
    Profile(ProfileArgs),
}

#[derive(ValueEnum, Clone, Copy, Debug)]
enum LibraryTypeArg {
    /// Read 1 is the reverse complement of the transcript, read 2 is on the
    /// transcript strand (dUTP, "fr-firststrand"; Salmon: ISR)
    Reverse,
    /// Read 1 is on the transcript strand, read 2 is its reverse complement
    /// ("fr-secondstrand"; Salmon: ISF)
    Forward,
    /// Either orientation with probability 1/2 (Salmon: IU)
    Unstranded,
}

#[derive(ValueEnum, Clone, Copy, Debug)]
enum AdapterFillArg {
    /// G, the "no signal" call of two-colour instruments
    #[value(name = "polyG")]
    PolyG,
    /// A
    #[value(name = "polyA")]
    PolyA,
    /// Uniformly random bases
    #[value(name = "random")]
    Random,
}

#[derive(ValueEnum, Clone, Copy, Debug)]
enum AbundanceTypeArg {
    /// Relative molar abundance (TPM or any value proportional to it):
    /// expected fragments are proportional to abundance x effective length
    Tpm,
    /// Expected fragments (read counts): expected fragments are proportional
    /// to the value itself
    Counts,
}

/// Parses a count such as `2000000`, `2_000_000`, `500k`, `2M` or `1.5M`.
fn parse_count(text: &str) -> std::result::Result<u64, String> {
    let cleaned: String = text.chars().filter(|c| *c != '_').collect();
    let (number, factor) = match cleaned.chars().last() {
        Some('k' | 'K') => (&cleaned[..cleaned.len() - 1], 1e3),
        Some('m' | 'M') => (&cleaned[..cleaned.len() - 1], 1e6),
        Some('g' | 'G') => (&cleaned[..cleaned.len() - 1], 1e9),
        _ => (cleaned.as_str(), 1.0),
    };
    if factor == 1.0 {
        return number.parse::<u64>().map_err(|_| {
            format!("'{text}' is not a whole number (suffixes k, M and G are allowed)")
        });
    }
    let value: f64 = number
        .parse()
        .map_err(|_| format!("'{text}' is not a number (suffixes k, M and G are allowed)"))?;
    let scaled = value * factor;
    if !scaled.is_finite() || scaled < 0.0 || scaled.fract() != 0.0 || scaled > 1e18 {
        return Err(format!("'{text}' is not a whole, non-negative number"));
    }
    Ok(scaled as u64)
}

#[derive(Args, Debug)]
struct SimulateArgs {
    /// Transcript sequences in FASTA format, plain or gzip. The identifier of
    /// a transcript is its header up to the first white space or '|'
    #[arg(short = 't', long, value_name = "FASTA", help_heading = "Input")]
    transcripts: PathBuf,

    /// Abundance table: tab-separated text with a header line, plain or gzip
    /// (for example a Salmon quant.sf, an RSEM isoforms.results or a kallisto
    /// abundance.tsv). Transcripts that are not in the table get abundance 0;
    /// identifiers that are not in the FASTA file are an error
    #[arg(short = 'a', long, value_name = "TSV", help_heading = "Input")]
    abundance: PathBuf,

    /// Name of the column with the transcript identifiers [default: the first
    /// of transcript_id, Name, target_id that exists]
    #[arg(long, value_name = "NAME", help_heading = "Input")]
    id_column: Option<String>,

    /// Name of the column with the abundances [default: tpm, in any case]
    #[arg(long, value_name = "NAME", help_heading = "Input")]
    abundance_column: Option<String>,

    /// What the abundances are
    #[arg(
        long,
        value_enum,
        default_value = "tpm",
        value_name = "TYPE",
        help_heading = "Input"
    )]
    abundance_type: AbundanceTypeArg,

    /// Prefix of the output files (may contain directories, which must exist)
    #[arg(short = 'o', long, value_name = "PREFIX", help_heading = "Output")]
    output_prefix: String,

    /// Write plain FASTQ files (<PREFIX>_1.fastq, <PREFIX>_2.fastq) instead
    /// of gzip files
    #[arg(long, help_heading = "Output")]
    uncompressed: bool,

    /// gzip compression level, 1 (fastest) to 9 (smallest files)
    #[arg(long, default_value_t = 4, value_name = "LEVEL",
          value_parser = clap::value_parser!(u32).range(1..=9), help_heading = "Output")]
    compression_level: u32,

    /// Also write <PREFIX>.reads.tsv.gz with the origin of every read pair
    /// (transcript, 1-based start, fragment length, strand of read 1)
    #[arg(long, help_heading = "Output")]
    per_read_truth: bool,

    /// Number of fragments (read pairs) to simulate; suffixes k, M and G are
    /// accepted, e.g. 2M
    #[arg(short = 'n', long, default_value = "1000000", value_name = "N",
          value_parser = parse_count, help_heading = "Library")]
    num_fragments: u64,

    /// Strandedness of the library
    #[arg(
        short = 'l',
        long,
        value_enum,
        default_value = "reverse",
        value_name = "TYPE",
        help_heading = "Library"
    )]
    library_type: LibraryTypeArg,

    /// Length of both reads, in bases
    #[arg(
        short = 'r',
        long,
        default_value_t = 100,
        value_name = "BASES",
        help_heading = "Library"
    )]
    read_length: usize,

    /// Mean of the normal fragment-length distribution before truncation, in
    /// bases
    #[arg(
        long,
        default_value_t = 180.0,
        value_name = "BASES",
        help_heading = "Library"
    )]
    fragment_mean: f64,

    /// Standard deviation of the fragment-length distribution before
    /// truncation, in bases (0 gives fragments of one length)
    #[arg(
        long,
        default_value_t = 55.0,
        value_name = "BASES",
        help_heading = "Library"
    )]
    fragment_sd: f64,

    /// Shortest fragment, in bases
    #[arg(
        long,
        default_value_t = 40,
        value_name = "BASES",
        help_heading = "Library"
    )]
    fragment_min: usize,

    /// Longest fragment, in bases
    #[arg(
        long,
        default_value_t = 600,
        value_name = "BASES",
        help_heading = "Library"
    )]
    fragment_max: usize,

    /// Sequence that read 1 runs into when the fragment is shorter than the
    /// read; {i7} and {i5} are replaced by the index sequences. The default
    /// is the Illumina TruSeq index adapter
    #[arg(long, default_value = DEFAULT_ADAPTER1, value_name = "SEQUENCE", help_heading = "Adapters")]
    adapter1: String,

    /// Sequence that read 2 runs into when the fragment is shorter than the
    /// read; {i7} and {i5} are replaced by the index sequences. The default
    /// is the reverse complement of the Illumina TruSeq dual-index i5 adapter
    #[arg(long, default_value = DEFAULT_ADAPTER2, value_name = "SEQUENCE", help_heading = "Adapters")]
    adapter2: String,

    /// Index 1 (i7) sequence: replaces {i7} in the adapters and is written
    /// in the read headers. May be empty
    #[arg(long, default_value = DEFAULT_INDEX1, value_name = "SEQUENCE", help_heading = "Adapters")]
    index1: String,

    /// Index 2 (i5) sequence, as read 2 reads it: replaces {i5} in the
    /// adapters and is written in the read headers after a '+'. May be empty
    /// (single-index library)
    #[arg(long, default_value = DEFAULT_INDEX2, value_name = "SEQUENCE", help_heading = "Adapters")]
    index2: String,

    /// Bases of a read past the end of the adapter
    #[arg(
        long,
        value_enum,
        default_value = "polyG",
        ignore_case = true,
        value_name = "FILL",
        help_heading = "Adapters"
    )]
    adapter_fill: AdapterFillArg,

    /// Quality profile written by the `profile` command [default: built-in
    /// profile of a recent instrument with binned qualities]
    #[arg(
        short = 'q',
        long,
        value_name = "JSON",
        help_heading = "Qualities and errors"
    )]
    quality_profile: Option<PathBuf>,

    /// Probability that a base is called N (with quality 2), between 0 and 1
    #[arg(
        long,
        default_value_t = 0.0,
        value_name = "PROBABILITY",
        help_heading = "Qualities and errors"
    )]
    n_rate: f64,

    /// Instrument name in the read names
    #[arg(
        long,
        default_value = "A00001",
        value_name = "NAME",
        help_heading = "Read names"
    )]
    instrument: String,

    /// Run number in the read names
    #[arg(
        long,
        default_value_t = 1,
        value_name = "NUMBER",
        help_heading = "Read names"
    )]
    run_number: u32,

    /// Flow cell name in the read names
    #[arg(
        long,
        default_value = "H7SIMDSXX",
        value_name = "NAME",
        help_heading = "Read names"
    )]
    flowcell: String,

    /// Lane number in the read names
    #[arg(
        long,
        default_value_t = 1,
        value_name = "NUMBER",
        help_heading = "Read names"
    )]
    lane: u32,

    /// Seed of the random number generator
    #[arg(
        short = 's',
        long,
        default_value_t = 1,
        value_name = "SEED",
        help_heading = "Run"
    )]
    seed: u64,

    /// Number of threads [default: all available cores]. The output does not
    /// depend on it
    #[arg(short = 'p', long, value_name = "THREADS", help_heading = "Run")]
    threads: Option<usize>,

    /// Do not print progress messages
    #[arg(long, help_heading = "Run")]
    quiet: bool,
}

#[derive(Args, Debug)]
struct ProfileArgs {
    /// FASTQ file of read 1 (or of single-end reads), plain or gzip
    #[arg(short = '1', long, value_name = "FASTQ")]
    read1: PathBuf,

    /// FASTQ file of read 2, plain or gzip. Without it the profile describes
    /// one read and `simulate` uses it for both reads
    #[arg(short = '2', long, value_name = "FASTQ")]
    read2: Option<PathBuf>,

    /// Profile file to write (JSON)
    #[arg(short = 'o', long, value_name = "JSON")]
    output: PathBuf,

    /// Number of reads to use from the start of each file; 0 uses all reads
    #[arg(long, default_value = "1000000", value_name = "N", value_parser = parse_count)]
    max_reads: u64,

    /// Treat cycles as independent: record only the per-cycle distributions,
    /// not the dependence on the quality of the previous cycle
    #[arg(long)]
    independent: bool,

    /// Do not print progress messages
    #[arg(long)]
    quiet: bool,
}

fn run_simulate(args: SimulateArgs) -> Result<()> {
    let threads = match args.threads {
        Some(0) => return Err(Error("--threads must be at least 1".into())),
        Some(n) => n,
        None => std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1),
    };
    let config = Config {
        transcripts: args.transcripts,
        abundance: args.abundance,
        id_column: args.id_column,
        abundance_column: args.abundance_column,
        abundance_type: match args.abundance_type {
            AbundanceTypeArg::Tpm => AbundanceType::Tpm,
            AbundanceTypeArg::Counts => AbundanceType::Counts,
        },
        output_prefix: args.output_prefix,
        quality_profile: args.quality_profile,
        reads: ReadParameters {
            num_fragments: args.num_fragments,
            seed: args.seed,
            library_type: match args.library_type {
                LibraryTypeArg::Reverse => LibraryType::Reverse,
                LibraryTypeArg::Forward => LibraryType::Forward,
                LibraryTypeArg::Unstranded => LibraryType::Unstranded,
            },
            read_length: args.read_length,
            fragment_mean: args.fragment_mean,
            fragment_sd: args.fragment_sd,
            fragment_min: args.fragment_min,
            fragment_max: args.fragment_max,
            adapter1: args.adapter1,
            adapter2: args.adapter2,
            index1: args.index1,
            index2: args.index2,
            adapter_fill: match args.adapter_fill {
                AdapterFillArg::PolyG => AdapterFill::PolyG,
                AdapterFillArg::PolyA => AdapterFill::PolyA,
                AdapterFillArg::Random => AdapterFill::Random,
            },
            n_rate: args.n_rate,
            instrument: args.instrument,
            run_number: args.run_number,
            flowcell: args.flowcell,
            lane: args.lane,
            per_read_truth: args.per_read_truth,
        },
        threads,
        compress: !args.uncompressed,
        compression_level: args.compression_level,
        quiet: args.quiet,
    };
    simulate::run(&config).map(|_| ())
}

fn run_profile(args: ProfileArgs) -> Result<()> {
    let log = |msg: &str| {
        if !args.quiet {
            eprintln!("[rnaseq-readsim] {msg}");
        }
    };
    let mut reads: Vec<ReadProfile> = Vec::new();
    let mut sources: Vec<String> = Vec::new();
    for path in std::iter::once(&args.read1).chain(args.read2.iter()) {
        let (profile, n) = learn_from_fastq(path, args.max_reads, !args.independent)?;
        log(&format!(
            "{}: {} reads, {} cycles, {} distinct quality values ({} to {})",
            path.display(),
            n,
            profile.cycles.len(),
            profile.qualities.len(),
            profile.qualities.first().copied().unwrap_or(0),
            profile.qualities.last().copied().unwrap_or(0),
        ));
        let file_name = path
            .file_name()
            .map(|f| f.to_string_lossy().into_owned())
            .unwrap_or_default();
        sources.push(format!("{file_name} ({n} reads)"));
        reads.push(profile);
    }
    let profile = Profile {
        format: PROFILE_FORMAT.to_owned(),
        version: PROFILE_VERSION,
        description: format!(
            "learned by rnaseq-readsim {} from {}; {}",
            env!("CARGO_PKG_VERSION"),
            sources.join(" and "),
            if args.independent {
                "independent cycles"
            } else {
                "first-order dependence on the previous cycle"
            }
        ),
        reads,
    };
    profile.validate().map_err(Error)?;
    std::fs::write(&args.output, profile.to_json())
        .context(|| format!("cannot write {}", args.output.display()))?;
    log(&format!("wrote {}", args.output.display()));
    Ok(())
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let result = match cli.command {
        Command::Simulate(args) => run_simulate(*args),
        Command::Profile(args) => run_profile(args),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::from(1)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_with_suffixes() {
        assert_eq!(parse_count("1000"), Ok(1000));
        assert_eq!(parse_count("2_000_000"), Ok(2_000_000));
        assert_eq!(parse_count("500k"), Ok(500_000));
        assert_eq!(parse_count("2M"), Ok(2_000_000));
        assert_eq!(parse_count("1.5m"), Ok(1_500_000));
        assert_eq!(parse_count("1G"), Ok(1_000_000_000));
        assert!(parse_count("1.5").is_err());
        assert!(parse_count("-3").is_err());
        assert!(parse_count("abc").is_err());
        assert!(parse_count("1.0001k").is_err());
        assert!(parse_count("").is_err());
    }

    #[test]
    fn command_line_is_consistent() {
        use clap::CommandFactory;
        Cli::command().debug_assert();
    }
}
