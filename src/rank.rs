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


#[cfg(test)]
mod test {
    use ndarray::array;

    use super::*;

    #[test]
    fn test_hrr() {
        assert_eq!(hrr(0, 1), 1);
        assert_eq!(hrr(1, 1), 1);
        assert_eq!(hrr(5, 1), 5);
    }

    #[test]
    fn test_mr() {
        assert_eq!(mr(1., 2.), (1.0 * 2.0).sqrt());
    }

    #[test]
    fn test_construct_rank_matrix_1() {
        let arr2 = array![
            [1.0, 0.9, 0.3],
            [0.9, 1.0, 0.5],
            [0.3, 0.5, 1.0]
        ];

        let rank: Array2<usize> = array![
            [0, 1, 2],
            [1, 0, 2],
            [2, 1, 0]
        ];

        assert_eq!(construct_rank_matrix(&arr2, 3).unwrap(), rank);
    }
}