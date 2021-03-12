use anyhow::Result;
use ndarray::Array2;
use ndarray_stats::*;
use std::path::PathBuf;

use crate::graph::Graph;
use crate::io;
use crate::rank;
use crate::Rank;

pub fn parse_args(
    input: &PathBuf,
    output: Option<&PathBuf>,
    method: Option<&Rank>,
    rank_cutoff: Option<&usize>,
    pcc_cutoff: Option<&f64>,
) -> Result<()> {
    info!("--- start read {}  ---", input.as_path().to_str().unwrap());

    // read csv and make ndarray::Array2
    let mut index: Vec<String> = vec![];

    let arr = io::read_exp_csv(input, &mut index)?;

    // calc correlation
    let corr = arr.pearson_correlation()?;
    debug!("{:?}", corr.shape());
    debug!("corr_matrix: \n{:?}", corr);
    debug!(
        "corr[0, 4]: {} {} : {:?}",
        &index[0],
        &index[4],
        corr[[0, 4]]
    );

    // calc rank matrix
    info!("calculate rank matrix...");
    let array_size = index.len();
    let rank_arr: Array2<usize> = rank::construct_rank_matrix(&corr, array_size)?;
    // construct hrr based network
    debug!("rank_matrix: \n{:?}", rank_arr);
    info!("construct rank based network...");

    let method = method.unwrap_or(&Rank::HRR);
    match method {
        Rank::HRR => {
            info!("Method: HRR");
            let mut g: Graph<usize> = Graph::new(&index);
            g.construct_hrr_network(corr, rank_arr, rank_cutoff, pcc_cutoff);
            let default_path = PathBuf::from("hrr_based_network.csv");
            let out_path = output.unwrap_or(&default_path);
            io::graph_to_csv(out_path.clone(), g)?;
        }
        Rank::MR => {
            info!("Method: MR");
            let mut g: Graph<f64> = Graph::new(&index);
            let rank_cutoff = match rank_cutoff {
                Some(rank_cutoff) => {
                    let rank_cutoff_f64 = *rank_cutoff as f64;
                    Some(rank_cutoff_f64)
                }
                None => None,
            };
            g.construct_mr_network(corr, rank_arr, rank_cutoff.as_ref(), pcc_cutoff);
            let default_path = PathBuf::from("mr_based_network.csv");
            let out_path = output.unwrap_or(&default_path);
            io::graph_to_csv(out_path.clone(), g)?;
        }
    }

    info!("Finish!");

    Ok(())
}
