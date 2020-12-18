use ndarray::{Array2};

use crate::io;
use crate::rank;

#[derive(Debug, Clone, Copy)]
pub struct Edge<T> 
    where T: Clone + Copy + PartialEq + PartialOrd
{
    node_1: usize,
    node_2: usize,
    corr: f64,
    rank: T,
}

impl<T> Edge<T> 
    where T: Clone + Copy + PartialEq + PartialOrd + ToString
{
    pub fn new(node_1: usize, node_2: usize, corr: f64, rank: T) -> Self {
        Self {
            node_1,
            node_2,
            corr,
            rank
        }
    }

    pub fn node_names(&self, nodes: &Vec<String>) -> (String, String) {
        (nodes[self.node_1].clone(), nodes[self.node_2].clone())
    }

    pub fn to_record(&self, nodes: &Vec<String>) -> io::CsvRecord {
        let (node_1_name, node_2_name) = self.node_names(nodes);
        io::CsvRecord::new(node_1_name, node_2_name, self.corr, self.rank.to_string())
    }
}

#[derive(Debug, Clone)]
struct Graph<T> 
    where T: Clone + Copy + PartialEq + PartialOrd
{
    edges: Vec<Edge<T>>,
    nodes: Vec<String>,
}

impl<T> Graph<T>
    where T: Clone + Copy + PartialEq + PartialOrd
{
    fn new(nodes: Vec<String>) -> Self {
        Self {
            edges: Vec::new(),
            nodes: nodes
        }
    }

    fn push(&mut self, edge: Edge<T>) {
        self.edges.push(edge)
    }

    fn size(&self) -> usize {
        self.nodes.len()
    }

    fn construct_hrr_network(
        &mut self,
        corr: Array2<f64>,
        rank: Array2<usize>,
        hrr_cutoff: usize,
        pcc_cutoff: Option<f64>,
    ) {
        for i in 0..self.size() {
            for j in i..self.size() {
                if let Some(pcc_cutoff) = pcc_cutoff {
                    if f64::abs(corr[[i, j]]) < pcc_cutoff { continue; }
                }

                let hrr = rank::hrr(rank[[i, j]], rank[[j, i]]);

                if hrr > hrr_cutoff { continue; }
                self.push(Edge::new(i, j, corr[[i, j]], hrr))
            }
        }
    }
}

