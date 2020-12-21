use std::{collections::HashMap, fmt::{Display, Debug}, str::FromStr, vec};
use ndarray::{Array2};

use crate::io;
use crate::rank;

#[derive(Debug, Clone)]
pub struct Node<T>
    where T: Clone + Copy
{
    node_name: String,
    edges: Vec<Edge<T>>,
}

impl<T> Node<T>
    where T: Clone + Copy
{
    pub fn new(node_name: String) -> Self {
        Self {
            node_name: node_name,
            edges: Vec::new()
        }
    }

    pub fn push(&mut self, edge: Edge<T>) {
        self.edges.push(edge)
    }

    pub fn edges(&self) -> Vec<Edge<T>> {
        self.edges.clone()
    }
}

impl<T> ToString for Node<T>
    where T: Clone + Copy
{
    fn to_string(&self) -> String {
        self.node_name.clone()
    }
}


#[derive(Debug, Clone, Copy)]
pub struct Edge<T> 
    where T: Clone + Copy
{
    node_1: usize,
    node_2: usize,
    corr: f64,
    rank: T,
}

impl<T> Edge<T> 
    where T: Clone + Copy
{
    pub fn new(node_1: usize, node_2: usize, corr: f64, rank: T) -> Self {
        Self {
            node_1,
            node_2,
            corr,
            rank
        }
    }

    pub fn query(&self) -> usize {
        self.node_1
    }

    pub fn target(&self) -> usize {
        self.node_2
    }

    pub fn node_names<S: ToString>(&self, nodes: &Vec<S>) -> (String, String) {
        (nodes[self.node_1].to_string(), nodes[self.node_2].to_string())
    }
}

impl<T> Edge<T>
    where T: Clone + Copy + ToString
{
    pub fn to_record<S: ToString>(&self, nodes: &Vec<S>) -> io::CsvRecord {
        let (node_1_name, node_2_name) = self.node_names(nodes);
        io::CsvRecord::new(node_1_name, node_2_name, self.corr, self.rank.to_string())
    }
}

#[derive(Debug, Clone)]
pub struct Graph<T> 
    where T: Clone + Copy
{
    // edges: Vec<Edge<T>>,
    nodes: Vec<Node<T>>,
}

impl<T> Graph<T>
    where T: Clone + Copy
{
    pub fn new(nodes: &Vec<String>) -> Self {
        Self {
            // edges: Vec::new(),
            nodes: nodes.iter().map(|x| Node::new(x.clone())).collect()
        }
    }

    fn push(&mut self, edge: Edge<T>) {
        let query = edge.query();
        self.nodes[query].push(edge);
    }

    fn size(&self) -> usize {
        self.nodes.len()
    }

    pub fn nodes(&self) -> &Vec<Node<T>> {
        &self.nodes
    }

    pub fn edges(&self) -> Vec<Edge<T>> {
        let mut edges = vec![];

        for n in self.nodes() {
            edges.extend(n.edges());
        }

        return edges
    }
}

// impl<T> Graph<T>
//     where T: Clone + Copy + FromStr + Debug + Display
// {
//     pub fn from_records(records: &Vec<io::CsvRecord>) -> Self
//         where <T as FromStr>::Err: std::fmt::Debug
//     {
//         let mut nodes: Vec<String> = vec![];
//         let mut edges: Vec<Edge<T>> = vec![];

//         let mut cnt = 0;
//         let mut map: HashMap<String, usize> = HashMap::new();

//         for r in records.iter() {
//             let (gene_1, gene_2) = r.genes();
//             let corr = r.corr();
//             let rank = r.rank::<T>();

//             let node_1 = get_node_index(&gene_1, &mut nodes, &mut map, &mut cnt);
//             let node_2 = get_node_index(&gene_2, &mut nodes, &mut map, &mut cnt);
//             edges.push(Edge::new(node_1, node_2, corr, rank));
//         }

//         Self {
//             edges,
//             nodes
//         }
//     }
// }


impl<T> Graph<T>
    where T: Clone + Copy + Ord + Display
{
    pub fn construct_hrr_network(
        &mut self,
        corr: Array2<f64>,
        rank: Array2<T>,
        rank_cutoff: Option<&T>,
        pcc_cutoff: Option<&f64>,
    ) {
        for i in 0..self.size() {
            for j in i..self.size() {
                if i == j { continue; }
                if let Some(pcc_cutoff) = pcc_cutoff {
                    if f64::abs(corr[[i, j]]) < *pcc_cutoff { continue; }
                }

                let hrr = rank::hrr(rank[[i, j]], rank[[j, i]]);

                if let Some(rank_cutoff) = rank_cutoff {
                    if hrr > *rank_cutoff { continue; }
                }
                self.push(Edge::new(i, j, corr[[i, j]], hrr))
            }
        }
    }
}

impl<T> Graph<T>
    where T: num_traits::Float + Display
{
    pub fn construct_mr_network<R: num_traits::NumCast + Copy>(
        &mut self,
        corr: Array2<f64>,
        rank: Array2<R>,
        rank_cutoff: Option<&T>,
        pcc_cutoff: Option<&f64>,
    ) {
        for i in 0..self.size() {
            for j in i..self.size() {
                if i == j { continue; }
                // construct network with rank::mr
                if let Some(pcc_cutoff) = pcc_cutoff {
                    if f64::abs(corr[[i, j]]) < *pcc_cutoff { continue; }
                }
                let mr = rank::mr(T::from(rank[[i, j]]).unwrap(), T::from(rank[[j, i]]).unwrap());
                
                if let Some(rank_cutoff) = rank_cutoff {
                    if mr > *rank_cutoff { continue; }
                }
                
                self.push(Edge::new(i, j, corr[[i, j]], mr));
            }
        }
    }
}

fn get_node_index(
    node_name: &String,
    nodes: &mut Vec<String>, 
    map: &mut HashMap<String, usize>,
    cnt: &mut usize
) -> usize {
    if map.contains_key(node_name) { return *map.get(node_name).unwrap() }

    let idx = *cnt;
    map.entry(node_name.to_owned()).or_insert(idx);
    nodes.push(node_name.to_owned());
    *cnt += 1;

    idx
}