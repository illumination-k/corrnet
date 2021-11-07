use std::{collections::HashMap, fmt::Debug};

use ndarray::{Array2};
use petgraph::{Graph, Undirected};

use crate::rank;

#[derive(Debug, Clone)]
pub struct Weight<T: Debug> {
    corr: f64,
    rank: T,
}

pub fn construct_hrr_network<T> (
    index: &[String],
    corr: &Array2<f64>,
    rank: &Array2<T>,
    rank_cutoff: Option<&T>,
    pcc_cutoff: Option<&f64>
) -> Graph<String, Weight<T>, Undirected> 
    where T: Ord + Copy +Debug
{
    let mut gr = Graph::new_undirected();
    let index_id2node_index = HashMap::new();
    for i in 0..index.len() {
        for j in i..index.len() {
            if i == j {
                continue;
            }

            if let Some(pcc_cutoff) = pcc_cutoff {
                if corr[[i, j]].abs() < *pcc_cutoff {
                    continue
                }
            }

            let hrr = rank::hrr(rank[[i, j]], rank[[j, i]]);

            if let Some(rank_cutoff) = rank_cutoff {
                if hrr > *rank_cutoff {
                    continue
                }
            }

            let w = Weight { corr: corr[[i, j]], rank: hrr };
            let node1 = gr.add_node(index[i].clone());
            let node2 = gr.add_node(index[j].clone());
            gr.add_edge(node1, node2, w);
        }
    }

    gr
}

#[cfg(test)]
mod test_construct_network {
    use super::*;
    use ndarray::array;
    use ndarray_stats::CorrelationExt;
    use crate::rank;

    #[test]
    fn test_hrr_network() {
        let index = ["gene_1", "gene_2", "gene_3"].map(|x| x.to_string());
        let arr2 = array![[1.0, 0.9, 0.3], [0.9, 1.0, 0.5], [0.3, 0.5, 1.0]];
        let corr = arr2.pearson_correlation().unwrap();
        let rank = rank::construct_rank_matrix(&corr, 3).unwrap();

        dbg!(&corr, &rank);
        let gr = construct_hrr_network(&index, &corr, &rank, None, None);
        dbg!(gr);
    }
}