extern crate anyhow;
extern crate csv;
extern crate bio;
extern crate pretty_env_logger;
extern crate ndarray;
extern crate ndarray_stats;
extern crate num_traits;

#[macro_use]
extern crate log;

extern crate serde;

#[macro_use]
extern crate serde_derive;

use std::{env::set_var};
use std::path::PathBuf;

use structopt::{clap, StructOpt, clap::arg_enum};
use anyhow::Result;

#[cfg(test)]
#[macro_use]
extern crate maplit;

mod io;
mod graph;
mod codon;
mod similarity;
mod rank;
mod handlers;
mod math;

#[derive(Debug, StructOpt)]
#[structopt(name = "hrr_corrnet")]
#[structopt(long_version(option_env!("LONG_VERSION").unwrap_or(env!("CARGO_PKG_VERSION"))))]
#[structopt(setting(clap::AppSettings::ColoredHelp))]
pub struct Opt {
    #[structopt(long = "log", possible_values(&LogLevel::variants()))]
    pub log_level: Option<LogLevel>,
    #[structopt(subcommand)]
    pub subcommand: SubCommands,
}

arg_enum! {
    #[derive(Debug)]
    pub enum LogLevel {
        DEBUG,
        INFO,
        WARN,
        ERROR,
    }
}

arg_enum! {
    #[derive(Debug)]
    pub enum Rank {
        HRR,
        MR,
    }
}

#[derive(Debug, StructOpt)]
pub enum SubCommands {
    #[structopt(name = "construct", about="construct rank based network from gene expression matrix")]
    #[structopt(setting(clap::AppSettings::ColoredHelp))]
    Construct {
        #[structopt(short = "-i", long = "input")]
        input: PathBuf,
        #[structopt(short = "-m", long = "method", possible_values(&Rank::variants()))]
        method: Option<Rank>,
        #[structopt(short = "-o", long = "output")]
        output: Option<PathBuf>,
        #[structopt(long = "rank_cutoff")]
        rank_cutoff: Option<usize>,
        #[structopt(long = "pcc_cutoff")]
        pcc_cutoff: Option<f64>
    },
    #[structopt(name = "extract", about="extract")]
    #[structopt(setting(clap::AppSettings::ColoredHelp))]
    Extract {
        #[structopt(short = "-i", long = "input")]
        input: PathBuf,
        #[structopt(short = "-g", long = "gene_list")]
        gene_list: Option<PathBuf>,
        #[structopt(short = "-o", long = "output")]
        output: Option<PathBuf>,
        #[structopt(long = "rank_cutoff")]
        rank_cutoff: Option<f64>,
        #[structopt(long = "pcc_cutoff")]
        pcc_cutoff: Option<f64>
    },
    #[structopt(name = "codon_usage", about="extract")]
    #[structopt(setting(clap::AppSettings::ColoredHelp))]
    CodonUsage {
        #[structopt(short = "-i", long = "input_graph")]
        input_graph: PathBuf,
        #[structopt(short = "-f", long = "input_fasta")]
        input_fasta: PathBuf,
        #[structopt(short = "-p", long = "percent")]
        percent: f64
    },
}

fn main() -> Result<()> {
    let opt = Opt::from_args();

    match &opt.log_level {
        Some(log_level) => {
            match log_level {
                LogLevel::DEBUG => set_var("RUST_LOG", "debug"),
                LogLevel::INFO => set_var("RUST_LOG", "info"),
                LogLevel::WARN => set_var("RUST_LOG", "warn"),
                LogLevel::ERROR => set_var("RUST_LOG", "error")
            }
        },
        None => set_var("RUST_LOG", "warn")
    };
    pretty_env_logger::init_timed();

    match &opt.subcommand {
        SubCommands::Construct {
            input,
            output,
            method,
            rank_cutoff,
            pcc_cutoff
        } => {
            handlers::construct::parse_args(
                input, 
                output.as_ref(), 
                method.as_ref(),
                rank_cutoff.as_ref(),
                pcc_cutoff.as_ref(),
            )?;
        },
        SubCommands::Extract {
            input,
            gene_list,
            output,
            rank_cutoff,
            pcc_cutoff
        } => {
            handlers::extract::parse_args(
                input, 
                gene_list.as_ref(), 
                output.as_ref(), 
                rank_cutoff.as_ref(), 
                pcc_cutoff.as_ref()
            )?;
        },
        SubCommands::CodonUsage {
            input_graph,
            input_fasta,
            percent
        } => {
            handlers::codon_usage::parse_args(
                input_graph,
                input_fasta,
                percent
            )?;
        }
    }
    Ok(())
}