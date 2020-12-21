use std::collections::HashSet;
use std::hash::Hash;

pub fn cosmix<T: Hash + Eq>(
    list: &Vec<T>,
    ref_list: &Vec<T>,
    k: usize
) -> f64 {
    // k should smaller than list.len() or equal
    assert!(k <= list.len());

    // list.len() == ref_list.len()
    assert!(list.len() == ref_list.len());

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
    use super::cosmix;

    #[test]
    fn test_cosmix() {
        let l = vec![0, 1, 2, 5, 6];
        let rl = vec![1, 2, 3, 4, 6];

        assert_eq!(cosmix(&l, &rl, 2), 1.0/3.0);
        assert_eq!(cosmix(&l, &rl, 3), 0.5);
        assert_eq!(cosmix(&l, &rl, 4), 5.0/10.0);
        assert_eq!(cosmix(&l, &rl, 5), 8.0/15.0);
    }
}