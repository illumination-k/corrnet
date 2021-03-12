use anyhow::Result;
use ndarray::{Array2, ArrayBase};
use num_traits::Float;
use ordered_float::OrderedFloat;
use superslice::*;

pub fn hrr<T: Ord>(a: T, b: T) -> T {
    std::cmp::max(a, b)
}

pub fn mr<T: Float>(a: T, b: T) -> T {
    (a * b).sqrt()
}

fn _rank(v: f64, sorted_vec: &Vec<OrderedFloat<f64>>) -> usize {
    let idx = sorted_vec.lower_bound(&OrderedFloat::from(v));
    sorted_vec.len() - idx
}

pub fn construct_rank_matrix(corr: &Array2<f64>, size: usize) -> Result<Array2<usize>> {
    let mut rank_vec = vec![];

    for row in corr.outer_iter() {
        let mut sorted_vec: Vec<OrderedFloat<f64>> = row
            .to_vec()
            .into_iter()
            .map(|x| OrderedFloat::from(x))
            .collect();
        sorted_vec.sort();
        for vv in row.to_vec().iter() {
            let rank = _rank(*vv, &sorted_vec);
            rank_vec.push(rank)
        }
    }

    Ok(ArrayBase::from_shape_vec((size, size), rank_vec)?)
}

pub fn get_index_sorted_by_rank(
    rank_matrix: &Array2<usize>,
    i: usize,
    index: &Vec<String>,
) -> Vec<String> {
    let mut rank_vec: Vec<String> = vec!["".to_string(); index.len() - 1];

    for j in 0..index.len() {
        let rank = rank_matrix[[i, j]];
        if rank == 0 {
            continue;
        }
        rank_vec[rank - 1] = index[j].clone();
    }

    rank_vec
}

#[cfg(test)]
mod test {
    use ndarray::array;

    use super::*;
    use crate::io;
    use ndarray_stats::*;

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
    fn test_rank_1() {
        let mut v: Vec<OrderedFloat<f64>> = [
            0.9999999999999999,
            -0.5788092174903143,
            0.23258203042199635,
            -0.6599591299048517,
            0.4039458299425666,
            -0.470197241576486,
            -0.7643118816809261,
            -0.392552428120889,
            -0.5913422295633524,
            -0.26489680851317693,
        ]
        .iter()
        .map(|x| OrderedFloat::from(*x))
        .collect();
        v.sort();
        assert_eq!(_rank(0.23258203042199635, &v), 2);
        assert_eq!(_rank(0.4039458299425666, &v), 1);
    }

    #[test]
    fn test_construct_rank_matrix_1() {
        let arr2 = array![[1.0, 0.9, 0.3], [0.9, 1.0, 0.5], [0.3, 0.5, 1.0]];

        let rank: Array2<usize> = array![[0, 1, 2], [1, 0, 2], [2, 1, 0]];

        assert_eq!(construct_rank_matrix(&arr2, 3).unwrap(), rank);
    }

    #[test]
    fn test_construct_rank_matrix_2() {
        let mut index: Vec<String> = vec![];
        let arr = io::read_exp_csv("test/small_test.csv", &mut index).unwrap();
        let corr = arr.pearson_correlation().unwrap();
        let rank = construct_rank_matrix(&corr, index.len()).unwrap();
        let expected: Array2<usize> = array![
            [0, 6, 2, 8, 1, 5, 9, 4, 7, 3],
            [9, 0, 4, 6, 8, 7, 5, 2, 1, 3],
            [3, 1, 0, 8, 6, 9, 7, 4, 2, 5],
            [8, 6, 9, 0, 3, 4, 2, 1, 5, 7],
            [1, 9, 6, 2, 0, 4, 5, 3, 8, 7],
            [8, 5, 9, 2, 4, 0, 1, 7, 6, 3],
            [9, 4, 8, 2, 7, 1, 0, 6, 5, 3],
            [9, 2, 7, 3, 5, 8, 6, 0, 1, 4],
            [9, 1, 4, 5, 8, 7, 6, 2, 0, 3],
            [8, 2, 7, 5, 9, 6, 1, 4, 3, 0]
        ];

        assert_eq!(rank, expected)
    }

    #[test]
    fn test_get_index_sorted_by_rank_1() {
        let rank: Array2<usize> = array![[0, 1, 2], [1, 0, 2], [2, 1, 0]];

        let index: Vec<String> = ["gene_1", "gene_2", "gene_3"]
            .iter()
            .map(|x| x.to_string())
            .collect();

        assert_eq!(
            get_index_sorted_by_rank(&rank, 0, &index),
            ["gene_2", "gene_3"]
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<String>>()
        );

        assert_eq!(
            get_index_sorted_by_rank(&rank, 1, &index),
            ["gene_1", "gene_3"]
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<String>>()
        );
        assert_eq!(
            get_index_sorted_by_rank(&rank, 2, &index),
            ["gene_2", "gene_1"]
                .iter()
                .map(|x| x.to_string())
                .collect::<Vec<String>>()
        );
    }
}
