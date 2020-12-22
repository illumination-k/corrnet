use std::collections::HashSet;
use std::hash::Hash;

use ordered_float::OrderedFloat;

pub fn median (list: &Vec<f64>) -> f64 {
    assert!(list.len() > 0);
    if list.len() == 1 { return list[0] }
    let mut v: Vec<OrderedFloat<f64>> = list.iter()
                                    .map(|x| OrderedFloat::from(*x))
                                    .collect();
    v.sort();

    if v.len() % 2 == 1 {
        v[v.len()/2].into()
    } else {
        ((v[v.len()/2 - 1] + v[v.len()/2]) / OrderedFloat::from(2.0)).into()
    }
}

pub fn cosmix<T: Hash + Eq>(
    list: &Vec<T>,
    ref_list: &Vec<T>,
    k: usize
) -> f64 {
    // k should smaller than list.len() or equal
    assert!(k <= list.len());
    // k should also smaller than ref_list.len() or equal
    assert!(k <= ref_list.len());

    let denominator: f64 = (1..=k).sum::<usize>() as f64;
    let mut numerator: f64 = 0.;
    for x in 1..=k {
        let set: HashSet<_> = list.iter().take(x).collect();
        let ref_set: HashSet<_> = ref_list.iter().take(x).collect();
        numerator += set.intersection(&ref_set).into_iter().count() as f64;
    }

    numerator / denominator
}


#[cfg(test)]
mod test {
    use std::vec;

    use super::{cosmix, median};

    #[test]
    fn test_cosmix_1() {
        let l = vec![0, 1, 2, 5, 6];
        let rl = vec![1, 2, 3, 4, 6];

        assert_eq!(cosmix(&l, &rl, 2), 1.0/3.0);
        assert_eq!(cosmix(&l, &rl, 3), 0.5);
        assert_eq!(cosmix(&l, &rl, 4), 5.0/10.0);
        assert_eq!(cosmix(&l, &rl, 5), 8.0/15.0);
    }

    #[test]
    fn test_cosmix_2() {
        let l = vec![1, 2, 3, 4, 5];
        let rl = vec![1, 2, 3, 4, 5];

        assert_eq!(cosmix(&l, &rl, 2), 3./3.);
        assert_eq!(cosmix(&l, &rl, 3), 6./6.);
    }

    #[test]
    fn test_median_1() {
        let odd = vec![1., 1., 2., 4., 5., 8., 9., 10., 11.];
        assert_eq!(median(&odd), 5.);

        let even = vec![1., 1., 2., 4., 5., 8., 9., 10., 11., 14.];
        assert_eq!(median(&even), 6.5);

        let one = vec![1.];
        assert_eq!(median(&one), 1.);
    }
}