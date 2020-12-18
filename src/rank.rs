pub fn hrr<T: Ord>(a: T, b: T) -> T {
    std::cmp::max(a, b)
}

pub fn mr(a: usize, b: usize) -> f64 {
    ((a as f64) * (b as f64)).sqrt()
}