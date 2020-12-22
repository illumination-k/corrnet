use std::path::PathBuf;
use std::collections::HashSet;

use csv::{Reader, Writer};
use anyhow::Result;

use crate::io;


pub fn parse_args(
    input: &PathBuf,
    gene_list: Option<&PathBuf>,
    output: Option<&PathBuf>,
    rank_cutoff: Option<&f64>,
    pcc_cutoff: Option<&f64>,
) -> Result<()> {
    let mut rdr = Reader::from_path(input)?;

    let default_path = PathBuf::from("extracted_network.csv");
    let out_path = output.unwrap_or(&default_path);
    let mut wtr = Writer::from_path(out_path)?;

    let gene_set: Option<HashSet<String>> = match gene_list {
        Some(p) => Some(io::read_gene_list(p)?),
        None => None
    };

    for _r in rdr.deserialize() {
        let r: io::CsvRecord = _r?;

        // filter by gene ids
        let (gene_1, gene_2) = r.genes();
        
        if let Some(gene_set) = gene_set.as_ref() {
            if !(gene_set.contains(&gene_1) && gene_set.contains(&gene_2)) {
                continue;
            }
        }

        // filter by rank
        if let Some(rank_cutoff) = rank_cutoff {
            let rank = r.rank::<f64>();
            if *rank_cutoff < rank { continue; }
        }

        // filter by pcc
        if let Some(pcc_cutoff) = pcc_cutoff {
            let pcc = r.corr();
            if pcc < *pcc_cutoff { continue; }
        }

        wtr.serialize(r)?;
    }

    wtr.flush()?;

    Ok(())
}