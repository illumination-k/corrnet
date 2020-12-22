use std::{collections::HashMap, path::PathBuf};
use anyhow::Result;
use csv::Reader;
use ordered_float::OrderedFloat;

use crate::io;
use crate::codon;
use crate::similarity;
use crate::rank;

pub fn parse_args(
    input_graph: &PathBuf,
    input_fasta: &PathBuf,
    percent: &f64,
) -> Result<()> {
    let (index, seqs) = io::read_fasta(input_fasta)?;

    // codon arr: convert to rank
    let codon_rank = codon::make_codon_rank(&seqs)?;


    let k = (index.len() as f64 * percent) as usize;

    // graph: sort by rank
    let mut rdr = Reader::from_path(input_graph)?;
    let mut map: HashMap<String, Vec<(String, OrderedFloat<f64>)>> = HashMap::new();

    for _r in rdr.deserialize() {
        let r: io::CsvRecord = _r?;
        let (gene_1, gene_2) = r.genes();
        let rank = OrderedFloat::from(r.rank::<f64>());
        map.entry(gene_1.clone())
            .or_insert(Vec::new())
            .push((gene_2.clone(), rank));
        map.entry(gene_2)
            .or_insert(Vec::new())
            .push((gene_1, rank));
    }

    // calc cosmix values
    let mut cosmix_values: Vec<f64> = vec![];

    for i in 0..index.len() {
        // let gene_id = index[i].clone();
        let corr_ranked_vec: Vec<String> = sort_corr_by_rank(&mut map, &index[i]);
        let codon_ranked_vec: Vec<String> = rank::get_index_sorted_by_rank(&codon_rank, i, &index);

        cosmix_values.push(
            similarity::cosmix(
                &corr_ranked_vec, 
                &codon_ranked_vec,
                k
            ) 
        );
    }

    // print median of comix values
    println!("Codon Score: {}", similarity::median(&cosmix_values));

    Ok(())
}

fn sort_corr_by_rank(
    map: &mut HashMap<String, Vec<(String, OrderedFloat<f64>)>>,
    key: &String
) -> Vec<String> {
    let corr_ranked = map.get_mut(key).unwrap();
    
    corr_ranked.sort_by(|a, b| a.1.cmp(&b.1));

    corr_ranked.iter()
        .map(|x| x.0.to_owned())
        .collect()
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_sort_corr_by_rank() {
        let mut map: HashMap<String, Vec<(String, OrderedFloat<f64>)>> = HashMap::new();
        let gene_1 = "gene_1".to_string();
        let gene_1_vec:  Vec<(String, OrderedFloat<f64>)>= vec![
            ("gene_2".to_string(), OrderedFloat::from(1.0)),
            ("gene_3".to_string(), OrderedFloat::from(4.0)),
            ("gene_4".to_string(), OrderedFloat::from(3.0)),
            ("gene_5".to_string(), OrderedFloat::from(2.0))
        ];
        map.insert(gene_1.clone(), gene_1_vec);

        let corr_ranked_vec = sort_corr_by_rank(&mut map, &gene_1);
        assert_eq!(
            corr_ranked_vec,
            ["gene_2", "gene_5", "gene_4", "gene_3"].iter()
                .map(|x| x.to_string())
                .collect::<Vec<String>>()
        )
    }
}