//! rnaseq-readsim: a fast, reproducible simulator of Illumina paired-end
//! bulk RNA-seq reads.
//!
//! The command-line program is a thin layer over this library:
//! [`simulate::run`] performs a whole run from files to files, and
//! [`simulate::Simulator`] produces chunks of reads in memory.

pub mod abundance;
pub mod alias;
pub mod error;
pub mod fasta;
pub mod fragment;
pub mod names;
pub mod quality;
pub mod reads;
pub mod rng;
pub mod simulate;
