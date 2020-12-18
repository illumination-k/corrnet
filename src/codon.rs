use std::collections::BTreeMap;
use anyhow::Result;
use ndarray::{Array2, ArrayBase};
use ndarray_stats::*;

fn make_codon_map() -> BTreeMap<String, usize> {
    let nuc = vec!['A', 'G', 'C', 'T'];
    let stop_codons = vec!["TGA", "TAA", "TAG"];
    let mut map: BTreeMap<String, usize> = BTreeMap::new();
    
    for n1 in nuc.iter() {
        for n2 in nuc.iter() {
            for n3 in nuc.iter() {
                let codon = format!("{}{}{}", n1, n2, n3);
                if stop_codons.contains(&codon.as_ref()) { continue; }
                map.entry(codon).or_insert(0);
            }
        }
    }
    map
}

fn make_codon_vec(
    seq: &str,
) -> Vec<usize> {
        let mut codon_map = make_codon_map();
        
        let seq: Vec<char> = seq.to_uppercase().chars().collect();

        for i in (0..seq.len()-3).step_by(3) {
            let codon: String = seq[i..i+3].into_iter().collect();
            match codon_map.get_mut(&codon) {
                Some(i) => *i += 1,
                None => continue
            }
        }

        codon_map.into_iter().map(|(_, x)| x).collect()
}

fn make_codon_arr(
    seqs: Vec<String>
) -> Result<Array2<usize>> {
    const codon_size: usize = 61;
    let seq_len = seqs.len();

    let mut vec: Vec<usize> = vec![];

    for seq in seqs.iter() {
        let codon_vec = make_codon_vec(seq);
        vec.extend(codon_vec);
    }

    Ok(ArrayBase::from_shape_vec((seq_len, codon_size), vec)?)
}

fn codon_corr(codon_arr: &Array2<usize>) -> Result<Array2<f64>> {
    Ok(codon_arr.mapv(|x| x as f64).pearson_correlation()?)
}