extern crate anyhow;
extern crate csv;
extern crate bio;
extern crate pretty_env_logger;
extern crate ndarray;
extern crate ndarray_stats;
extern crate chrono;

#[macro_use]
extern crate log;

use std::env::set_var;
use std::path::{Path, PathBuf};

use chrono::{Local};

use ndarray::{ArrayBase, Axis, Array1, Array2, aview0};
use ndarray_stats::*;

use ordered_float::OrderedFloat;
use superslice::*;

use structopt::{clap, StructOpt, clap::arg_enum};
use anyhow::Result;

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

    let mut rdr = csv::Reader::from_path(&opt.input)?;
    info!("--- start read {}: {} ---", &opt.input.as_path().to_str().unwrap(), Local::now());

    // read csv and make ndarray::Array2
    let mut index: Vec<String> = vec![];
    let mut vec: Vec<f64> = vec![];

    let mut shape = (0, 0);
    for _r in rdr.records() {
        let r = _r?;
        
        // skip index
        let exp_vec: Vec<f64> = r.into_iter()
            .skip(1)
            .map(|x| x.parse::<f64>().expect("cannot convert to f64"))
            .collect();
        
        // Whine std == 0, continue. Maybe use approx for abs_diff_eq
        let exp_arr: Array1<f64> = ArrayBase::from(exp_vec.clone());
        if exp_arr.std_axis(Axis(0), 1.0)  == aview0(&0.) { continue; }

        index.push(r[0].to_string());
        vec.extend(exp_vec);
        shape.1 = r.len()-1;
        shape.0 += 1;
        // debug!("{:?}", &r)
    }
    info!("shape: {:?}", shape);
    info!("caluclate pearson correlation... : {}", Local::now());
    let arr: Array2<f64> = ArrayBase::from_shape_vec(shape, vec)?;

    // calc correlation
    let corr = arr.pearson_correlation()?;
    debug!("{:?}", corr.shape());
    debug!("corr_matrix: \n{:?}", corr);
    debug!("corr[0, 4]: {} {} : {:?}", &index[0], &index[4], corr[[0, 4]]);

    // calc rank matrix
    info!("calculate rank matrix... : {}", Local::now());
    let mut rank_vec: Vec<usize> = vec![];
    for row in corr.outer_iter() {
        // sort and bitsect
        let mut sorted_vec: Vec<OrderedFloat<f64>> = row.to_vec().into_iter().map(|x| OrderedFloat::from(f64::abs(x))).collect();
        sorted_vec.sort();

        for vv in row.to_vec().iter() {
            let rank = sorted_vec.len() - sorted_vec.lower_bound(&OrderedFloat::from(*vv)) - 1;
            rank_vec.push(rank)
        }
    }

    let array_size = index.len();
    let rank_arr: Array2<usize>  = ArrayBase::from_shape_vec((array_size, array_size), rank_vec)?;

    // construct hrr based network
    info!("construct hrr based network... : {}", Local::now());
    let hrr_cutoff: usize = *(&opt.hrr_cutoff);
    let pcc_cutoff: Option<f64> = *(&opt.pcc_cutoff);
    let mut graph = Graph::new(index);

    for i in 0..array_size {
        for j in i..array_size {
            if i == j { continue; }
            if let Some(pcc_cutoff) = pcc_cutoff {
                if f64::abs(corr[[i, j]]) < pcc_cutoff { continue; }
            }
            // if std::cmp::min(rank_arr[[i, j]], rank_arr[[j, i]]) > hrr_cutoff { continue; }
            // highest reciprocal rank: max(rank(A, B), rank(B, A))
            let hrr = std::cmp::max(rank_arr[[i, j]], rank_arr[[j, i]]);

            if hrr > hrr_cutoff { continue; }
            // println!("query: {} target: {} corr: {} hrr: {}", &graph.index[i], &graph.index[j], corr[[i, j]], hrr);
            graph.push(Node::new(i, j, corr[[i, j]], hrr));
        }
    }

    // calc hcca clusters

    info!("Write csv...: {}", Local::now());
    graph_to_csv("ignore/graph.csv", graph)?;
    // calc codon usage
    Ok(())
}

fn graph_to_csv<P: AsRef<Path>>(outpath: P, graph: Graph) -> Result<()> {
    let mut wtr = csv::Writer::from_path(outpath.as_ref())?;
    wtr.write_record(&["query", "target", "corr", "hrr"])?;

    // use serializer ?
    for node in graph.nodes.iter() {
        wtr.write_record(&[node.query_name(&graph.index), node.target_name(&graph.index), node.corr.to_string(), node.hrr.to_string()])?;
    }

    wtr.flush()?;

    Ok(())
}

#[derive(Debug, Clone, Copy)]
struct Node {
    query: usize,
    target: usize,
    corr: f64,
    hrr: usize,
}

impl Node {
    fn new(query: usize, target: usize, corr: f64, hrr: usize) -> Self {
        Self {
            query,
            target,
            corr,
            hrr
        }
    }

    fn query_name(&self, index: &Vec<String>) -> String {
        index[self.query].clone()
    }

    fn target_name(&self, index: &Vec<String>) -> String {
        index[self.target].clone()
    }
}

#[derive(Debug, Clone)]
struct Graph {
    nodes: Vec<Node>,
    index: Vec<String>,
}

impl Graph {
    fn new(index: Vec<String>) -> Self {
        Self {
            nodes: Vec::new(),
            index: index
        }
    }

    fn push(&mut self, node: Node) {
        self.nodes.push(node)
    }
}