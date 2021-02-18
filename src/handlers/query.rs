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
    let mut wtr = Writer::from_path("test_query.csv")?;

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
            let byte_rec = io::CsvRecord::from_byte_records(&record);
            wtr.serialize(byte_rec)?;
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

    let mut edges = vec![];
    dfs(gene_id, 0, depth, &mut edges, &graph);

    for edge in edges.iter() {
        let rec = io::CsvRecord::from_tuple(edge);
        wtr.serialize(rec)?;
    }
    wtr.flush()?;
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
        let (corr, rank) = match map.get(k) {
            Some((corr, rank)) => { (corr, rank) },
            None => {continue;}
        };
        edges.push((query.clone(), k.clone(), *corr, *rank));
        dfs(k, depth + 1, depth_limit, edges, graph);
    }
}

#[cfg(test)]
pub mod test {
    use std::collections::HashMap;

    use super::dfs;

    #[test]
    fn test_dfs() {
        let edges = vec![
            ("gene_1", "gene_2", 1., 1.),
            ("gene_1", "gene_3", 1., 1.),
            ("gene_2", "gene_6", 1., 1.),
            ("gene_3", "gene_4", 1., 1.,),
            ("gene_4", "gene_5", 1., 1.,),
        ];

        let mut graph: HashMap<String, HashMap<String, (f64, f64)>> = HashMap::new();

        for edge in edges.iter() {
            graph.entry(edge.0.to_string()).or_insert(HashMap::new())
                .entry(edge.1.to_string()).or_insert((edge.2, edge.3));
        }

        let mut dfs_edges = vec![];
        dfs(&"gene_1".to_string(), 0, 2, &mut dfs_edges, &graph);

        let ans: Vec<(String, String, f64, f64)> = vec![
            ("gene_1", "gene_2", 1., 1.),
            ("gene_1", "gene_3", 1., 1.),
            ("gene_2", "gene_6", 1., 1.),
            ("gene_3", "gene_4", 1., 1.,)
        ].iter()
            .map(|(g1, g2, corr, rank)| (g1.to_string(), g2.to_string(), *corr, *rank))
            .collect();
        
        for edge in dfs_edges.iter() {
            assert!(ans.contains(edge))
        }
    }
}