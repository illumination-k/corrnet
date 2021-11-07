use std::{ffi::OsStr, fs::File, path::PathBuf};

use anyhow::Result;
use flate2::{Compression, GzBuilder};
use polars::prelude::*;

use crate::Rank;

fn read_graph_by_poloars(path: &PathBuf, rank: &Rank) -> Result<DataFrame> {
    let schema = Schema::new(vec![
        Field::new("gene_1", DataType::Utf8),
        Field::new("gene_2", DataType::Utf8),
        Field::new("corr", DataType::Float64),
        Field::new("rank", DataType::Float64),
    ]);

    let mut df = CsvReader::from_path(path)?
        .has_header(true)
        .with_schema(&schema)
        .finish()?;

    let col_name = match rank {
        Rank::HRR => "hrr_rank",
        Rank::MR => "mr_rank",
    };

    df.rename("rank", col_name)?;

    Ok(df)
}

fn add_relation(df: &mut DataFrame) -> Result<()> {
    let gene1chunk = df.column("gene_1")?.utf8()?;
    // Either ownership is required for add op with "-"
    let gene2chunk = df.column("gene_2")?.utf8()?.to_owned();
    let genechunk = gene1chunk + "-" + gene2chunk.to_owned();
    df.replace_or_add("relation", genechunk.into_series())?;
    Ok(())
}

fn filterout_by_rank_and_merge(
    df: &DataFrame,
    other: &DataFrame,
    priority: &Rank,
    max_rank: f64,
) -> Result<DataFrame> {
    let col_name = match priority {
        Rank::HRR => "hrr_rank",
        Rank::MR => "mr_rank",
    };

    let mask = df.column(col_name)?.f64()?.lt_eq(max_rank);
    let df = df.filter(&mask)?;

    let merged_df = df
        .inner_join(other, "relation", "relation")?
        .select(vec!["gene_1", "gene_2", "corr", "hrr_rank", "mr_rank"])?;
    Ok(merged_df)
}

pub fn parse_args(
    hrr_path: &PathBuf,
    mr_path: &PathBuf,
    out_path: &PathBuf,
    priority: &Rank,
    max_rank: &f64,
) -> Result<()> {
    info!(
        "\n hrr graph path: {:?}\n mr graph path: {:?}\n priority rank: {}, max rank: {}",
        hrr_path, mr_path, priority, max_rank
    );

    let mut hrr_df = read_graph_by_poloars(hrr_path, &Rank::HRR)?;
    let mut mr_df = read_graph_by_poloars(mr_path, &Rank::MR)?;

    add_relation(&mut hrr_df)?;
    add_relation(&mut mr_df)?;

    let merged_df = match priority {
        Rank::HRR => filterout_by_rank_and_merge(&hrr_df, &mr_df, priority, *max_rank),
        Rank::MR => filterout_by_rank_and_merge(&mr_df, &hrr_df, priority, *max_rank),
    }?;

    let ext = out_path.extension();

    let mut f = File::create(out_path)?;
    if ext == Some(OsStr::new("gz")) {
        let mut gz = GzBuilder::new().write(f, Compression::default());
        CsvWriter::new(&mut gz)
            .has_header(true)
            .with_delimiter(b',')
            .finish(&merged_df)?;
    } else {
        CsvWriter::new(&mut f)
            .has_header(true)
            .with_delimiter(b',')
            .finish(&merged_df)?;
    }

    Ok(())
}
