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
    let mut rdr = Reader::from_path("test")?;
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
                if record.corr() < pcc_cutoff { continue; }
            }

            if let Some(rank_cutoff) = rank_cutoff {
                if record.rank() > rank_cutoff { continue; }
            }

            wtr.serialize(record)?;
        }

        wtr.flush()?;
        return Ok(())
    }

    // 一端全部読み込むけど、pccとかrankとかがダメな奴は無視していいはず
    
    Ok(())
}