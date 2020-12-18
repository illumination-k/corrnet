#[derive(Debug, Clone, Copy)]
pub struct Edge {
    node_1: usize,
    node_2: usize,
    corr: f64,
    hrr: usize,
}

impl Edge {
    pub fn new(node_1: usize, node_2: usize, corr: f64, hrr: usize) -> Self {
        Self {
            node_1,
            node_2,
            corr,
            hrr
        }
    }
}

