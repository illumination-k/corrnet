use std::path::PathBuf;

use anyhow::Result;

<<<<<<< HEAD
#[allow(dead_code)]
pub fn parse_args(
    _in_graph: &PathBuf,
    _rank_cutoff: Option<&f64>,
=======

pub fn parse_args(
    in_graph: &PathBuf,
    rank_cutoff: Option<&f64>,
>>>>>>> origin/main
) -> Result<()> {
    Ok(())
}