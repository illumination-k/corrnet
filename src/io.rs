use std::path::Path;
use std::io::{BufRead, BufReader};
use std::collections::HashMap;

use chrono::{Local};
use flate2::read::MultiGzDecoder;

use csv::{Reader, Writer};
use anyhow::Result;
use ndarray::{Axis, ArrayBase, Array1, Array2, aview0};

use crate::graph;

pub fn read_exp_csv<P: AsRef<Path>>(
    input: P, 
    index: &mut Vec<String>, 
) -> Result<Array2<f64>> {
    let mut shape = (0, 0);
    let mut vec: Vec<f64> = vec![];
    let mut rdr = Reader::from_path(input)?;

    for _r in rdr.records() {
        let r = _r?;
        // skip index
        let exp_vec: Vec<f64> = r.into_iter()
            .skip(1)
            .map(|x| x.parse::<f64>().expect("cannot convert to f64"))
            .collect();
        
        // Whine std == 0, continue. Maybe use approx for abs_diff_eq
        let exp_arr: Array1<f64> = ArrayBase::from(exp_vec.clone());
        if exp_arr.std_axis(Axis(0), 1.0)  == aview0(&0.) { continue; }

        index.push(r[0].to_string());
        vec.extend(exp_vec);
        shape.1 = r.len()-1;
        shape.0 += 1;
        
    }
    info!("shape: {:?}", shape);
    info!("caluclate pearson correlation... : {}", Local::now());
    Ok(ArrayBase::from_shape_vec(shape, vec)?)
}

fn open_with_gz<P: AsRef<Path>>(p: P) -> Result<Box<dyn BufRead>> {
    let r = std::fs::File::open(p.as_ref())?;
    let ext = p.as_ref().extension();

    if ext == Some(std::ffi::OsStr::new("gz")) {
        let gz = MultiGzDecoder::new(r)?;
        let buf_reader = BufReader::new(gz);
        Ok(Box::new(buf_reader))
    } else {
        let buf_reader = BufReader::new(r);
        Ok(Box::new(buf_reader))
    }
}

pub fn read_fasta<P: AsRef<Path>>(fasta: P) -> Result<HashMap<String, String>> {
    let mut map: HashMap<String, String> = HashMap::new();
    let rdr = bio::io::fasta::Reader::new(open_with_gz(fasta)?);
    
    for _r in rdr.records() {
        let r = _r?;
        map.entry(r.id().to_string()).or_insert(String::from_utf8(r.seq().to_vec())?);
    }

    Ok(map)
}

#[derive(Debug, Serialize)]
pub struct CsvRecord {
    gene_1: String,
    gene_2: String,
    corr: f64,
    rank: String,
}

impl CsvRecord {
    pub fn new(gene_1: String, gene_2: String, corr: f64, rank: String) -> Self {
        Self {
            gene_1,
            gene_2,
            corr,
            rank
        }
    }
}