extern crate anyhow;
extern crate csv;
extern crate bio;
extern crate pretty_env_logger;
extern crate ndarray;
extern crate ndarray_stats;
extern crate chrono;
extern crate num_traits;

#[macro_use]
extern crate log;

extern crate serde;

#[macro_use]
extern crate serde_derive;

use std::env::set_var;
use std::path::PathBuf;

use chrono::{Local};

use ndarray::Array2;
use ndarray_stats::*;

use structopt::{clap, StructOpt, clap::arg_enum};
use anyhow::Result;

mod io;
mod graph;
mod codon;
mod similarity;
mod rank;

#[derive(Debug, StructOpt)]
#[structopt(name = "hrr_corrnet")]
#[structopt(long_version(option_env!("LONG_VERSION").unwrap_or(env!("CARGO_PKG_VERSION"))))]
#[structopt(setting(clap::AppSettings::ColoredHelp))]
pub struct Opt {
    #[structopt(short = "i", long = "input")]
    pub input: PathBuf,
    #[structopt(short = "t", long = "threads", default_value("1"), value_name = "INT")]
    pub threads: u8,
    #[structopt(long = "log", possible_values(&LogLevel::variants()))]
    pub log_level: Option<LogLevel>,
    #[structopt(long = "hrr_cutoff", default_value("30"), value_name = "INT")]
    pub hrr_cutoff: usize,
    #[structopt(long = "pcc_cutoff", value_name = "FLOAT")]
    pub pcc_cutoff: Option<f64>,
    #[structopt(short = "-o", long = "out_dir", default_value("hrr_results"))]
    pub outdir: PathBuf,
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
    pretty_env_logger::init();

    info!("--- start read {}: {} ---", &opt.input.as_path().to_str().unwrap(), Local::now());

    // read csv and make ndarray::Array2
    let mut index: Vec<String> = vec![];
    
    let arr = io::read_exp_csv(&opt.input, &mut index)?;

    // calc correlation
    let corr = arr.pearson_correlation()?;
    debug!("{:?}", corr.shape());
    debug!("corr_matrix: \n{:?}", corr);
    debug!("corr[0, 4]: {} {} : {:?}", &index[0], &index[4], corr[[0, 4]]);

    // calc rank matrix
    info!("calculate rank matrix... : {}", Local::now());
    let array_size = index.len();
    let rank_arr: Array2<usize> = rank::construct_rank_matrix(&corr, array_size)?;
    // construct hrr based network
    info!("construct hrr based network... : {}", Local::now());
    let hrr_cutoff: usize = *(&opt.hrr_cutoff);
    let pcc_cutoff: Option<f64> = *(&opt.pcc_cutoff);
    let mut graph: graph::Graph<usize> = graph::Graph::new(&index);
    graph.construct_hrr_network(corr, rank_arr, hrr_cutoff, pcc_cutoff);

    // calc hcca clusters

    info!("Write csv...: {}", Local::now());
    io::graph_to_csv("ignore/graph.csv", graph)?;
    // calc codon usage
    Ok(())
}