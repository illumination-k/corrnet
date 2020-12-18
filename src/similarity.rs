use std::collections::HashSet;
use std::hash::Hash;

pub fn cosmix<T: Hash + Eq>(
    list: Vec<T>,
    ref_list: Vec<T>,
    k: usize
) -> f64 {
    let denominator: f64 = (1..=k).sum::<usize>() as f64;
    let mut numerator: f64 = 0.;
    for x in 1..=k {
        let set: HashSet<_> = list.iter().take(x).collect();
        let ref_set: HashSet<_> = ref_list.iter().take(x).collect();
        numerator += set.intersection(&ref_set).into_iter().count() as f64;
    }

    numerator / denominator
}