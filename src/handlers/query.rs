use std::{collections::HashMap, path::PathBuf};

use anyhow::Result;
use csv::{Reader, Writer};

use crate::io;

pub fn parse_args(
    gene_id: &String,
    input_path: &PathBuf,
    depth: usize,
    pcc_cutoff: Option<&f64>,
    rank_cutoff: Option<&f64>,
) -> Result<()> {
    let gene_id_bytes = gene_id.as_bytes();
    let mut rdr = Reader::from_path(input_path)?;
    let mut wtr = Writer::from_path("test")?;

    let mut raw_record = csv::ByteRecord::new();
    let headers = rdr.byte_headers()?.clone();

    if depth == 1 {
        while rdr.read_byte_record(&mut raw_record)? {
            let record: io::ByteCsvRecord = raw_record.deserialize(Some(&headers))?;
            if record.gene_1_bytes() != gene_id_bytes && record.gene_2_bytes() != gene_id_bytes {
                continue;
            }

            if let Some(pcc_cutoff) = pcc_cutoff {
                if record.corr() < *pcc_cutoff { continue; }
            }

            if let Some(rank_cutoff) = rank_cutoff {
                if record.rank() > *rank_cutoff { continue; }
            }

            wtr.serialize(record)?;
        }

        wtr.flush()?;
        return Ok(())
    }

    // 一端全部読み込むけど、pccとかrankとかがダメな奴は無視していいはず
    // とりあえずhashmapのグラフでdfsか
    let mut graph: HashMap<String, HashMap<String, (f64, f64)>> = HashMap::new();
    
    while rdr.read_byte_record(&mut raw_record)? {
        let record: io::ByteCsvRecord = raw_record.deserialize(Some(&headers))?;


        if let Some(pcc_cutoff) = pcc_cutoff {
            if record.corr() < *pcc_cutoff { continue; }
        }

        if let Some(rank_cutoff) = rank_cutoff {
            if record.rank() > *rank_cutoff { continue; }
        }
        
        let (gene_1, gene_2) = record.genes_unchecked();

        graph.entry(gene_1.clone()).or_insert(HashMap::new());
        graph.get_mut(&gene_1)
            .unwrap()
            .entry(gene_2.clone())
            .or_insert((record.corr(), record.rank()));

        graph.entry(gene_2.clone()).or_insert(HashMap::new());
        graph.get_mut(&gene_2)
            .unwrap()
            .entry(gene_1)
            .or_insert((record.corr(), record.rank()));
    }

    Ok(())
}

fn dfs(
    query: &String,
    depth: usize,
    depth_limit: usize,
    edges: &mut Vec<(String, String, f64, f64)>,
    graph: &HashMap<String, HashMap<String, (f64, f64)>>
) {
    if depth == depth_limit { return; }
    let map = match graph.get(query) {
        Some(map) => map,
        None => return
    };

    for k in map.keys() {
        let (corr, rank) = map.get(query).unwrap();
        edges.push((query.clone(), k.clone(), *corr, *rank));
        dfs(k, depth + 1, depth_limit, edges, graph);
    }
}

#[cfg(test)]
pub mod test {
    use super::dfs;
}