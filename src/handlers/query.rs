use anyhow::Result;
use csv::{Reader, Writer};

use crate::io;

pub fn parse_args(
    gene_id: &str,
    depth: usize,
    pcc_cutoff: Option<f64>,
    rank_cutoff: Option<f64>,
) -> Result<()> {
    let gene_id_bytes = gene_id.as_bytes();

    let mut wtr = Writer::from_path("test")?;
    if depth == 1 {
        return Ok(())
    }


    
    Ok(())
}