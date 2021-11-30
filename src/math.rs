use ordered_float::OrderedFloat;

pub fn mean(list: &[f64]) -> f64 {
    list.iter().sum::<f64>() / list.len() as f64
}

pub fn var(list: &[f64], ddof: f64) -> f64 {
    let mean = mean(list);
    list.iter().map(|x| (x - mean).powi(2i32)).sum::<f64>() / (list.len() as f64 - ddof)
}

pub fn std(list: &[f64], ddof: f64) -> f64 {
    var(list, ddof).sqrt()
}

pub fn median(list: &[f64]) -> f64 {
    assert!(!list.is_empty());
    if list.len() == 1 {
        return list[0];
    }
    let mut v: Vec<OrderedFloat<f64>> = list.iter().map(|x| OrderedFloat::from(*x)).collect();
    v.sort();

    if v.len() % 2 == 1 {
        v[v.len() / 2].into()
    } else {
        ((v[v.len() / 2 - 1] + v[v.len() / 2]) / OrderedFloat::from(2.0)).into()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use approx::*;

    #[test]
    fn test_median_1() {
        let odd = vec![1., 1., 2., 4., 5., 8., 9., 10., 11.];
        assert_eq!(median(&odd), 5.);

        let even = vec![1., 1., 2., 4., 5., 8., 9., 10., 11., 14.];
        assert_eq!(median(&even), 6.5);

        let one = vec![1.];
        assert_eq!(median(&one), 1.);

        let rand10 = vec![
            0.32840955, 0.48140666, 0.1176708, 0.10189263, 0.53973073, 0.49730681, 0.42883597,
            0.86240549, 0.84503774, 0.22184689,
        ];
        assert_abs_diff_eq!(median(&rand10), 0.45512131499999997);

        let rand9 = vec![
            0.32840955, 0.48140666, 0.1176708, 0.10189263, 0.53973073, 0.49730681, 0.42883597,
            0.86240549, 0.84503774,
        ];
        assert_eq!(median(&rand9), 0.48140666);
    }

    #[test]
    fn test_std_1() {
        let rand5 = vec![0.30330361, 0.04612777, 0.41467306, 0.15042536, 0.01180612];
        assert_abs_diff_eq!(std(&rand5, 1.), 0.171188582970728);
        assert_abs_diff_eq!(std(&rand5, 0.), 0.1531157233977643);

        let rand10 = vec![
            0.32840955, 0.48140666, 0.1176708, 0.10189263, 0.53973073, 0.49730681, 0.42883597,
            0.86240549, 0.84503774, 0.22184689,
        ];
        assert_abs_diff_eq!(std(&rand10, 1.), 0.265779165304154);
        assert_abs_diff_eq!(std(&rand10, 0.), 0.2521402550938575);
    }
}
