pub fn hrr(a: usize, b: usize) -> usize {
    std::cmp::max(a, b)
}

pub fn mr(a: usize, b: usize) -> f64 {
    ((a as f64) * (b as f64)).sqrt()
}