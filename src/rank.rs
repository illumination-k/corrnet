use anyhow::Result;
use num_traits::Float;
use ndarray::{Array2, ArrayBase};
use ordered_float::OrderedFloat;
use superslice::*;

pub fn hrr<T: Ord>(a: T, b: T) -> T {
    std::cmp::max(a, b)
}

pub fn mr<T: Float>(a: T, b: T) -> T {
    (a * b).sqrt()
}


pub fn construct_rank_matrix(corr: &Array2<f64>, size: usize) -> Result<Array2<usize>> {
    let mut rank_vec = vec![];

    for row in corr.outer_iter() {
        let mut sorted_vec: Vec<OrderedFloat<f64>> = row.to_vec()
                                                            .into_iter()
                                                            .map(|x| OrderedFloat::from(f64::abs(x)))
                                                            .collect();
        sorted_vec.sort();

        for vv in row.to_vec().iter() {
            let rank = sorted_vec.len() - sorted_vec.lower_bound(&OrderedFloat::from(*vv)) - 1;
            rank_vec.push(rank)
        }
    }

    Ok(ArrayBase::from_shape_vec((size, size), rank_vec)?)
}