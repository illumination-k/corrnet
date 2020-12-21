use std::{collections::HashMap, path::PathBuf};
use anyhow::Result;
use csv::Reader;
use ordered_float::OrderedFloat;

use crate::io;
use crate::codon;
use crate::similarity;

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
        let gene_id = index[i].clone();
        let corr_ranked = map.get_mut(&gene_id).unwrap();
        corr_ranked.sort_by(|a, b| a.1.cmp(&b.1));
        let corr_ranked_vec: Vec<String> = corr_ranked.iter().map(|x| x.0.to_owned()).collect();
        let mut codon_ranked_vec: Vec<String> = vec!["".to_string(); index.len() - 1];

        for j in 0..index.len()-1 {
            let rank = codon_rank[[i, j]];
            if rank == 0 { continue; }
            codon_ranked_vec[rank-1] = index[j].clone();
        }

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