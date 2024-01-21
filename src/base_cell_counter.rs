use anyhow::Result;
use clap::{command, Parser};
use noodles::bed::{self, Record};
use rust_htslib::faidx;
use std::{cmp, fs::File, io::BufReader, time::Instant};

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    #[arg(long)]
    fasta: String,
    #[arg(long)]
    chrom: String,
    #[arg(long)]
    bin: usize,
    #[arg(long)]
    bed: Option<String>,
    #[arg(long)]
    bed_out: Option<String>,
}

pub fn main() -> Result<()> {
    let args = &Args::parse();

    let start = Instant::now();

    let _res = make_windows(args)?;

    let elapsed = start.elapsed();
    println!("Computation time: {} ms", elapsed.as_millis());

    Ok(())
}

fn make_windows(args: &Args) -> Result<Vec<Record<3>>> {
    // Pass 1. Get bed file and focus on the selected chromosome (if specified)
    let a = get_records(args)?;

    // Focus only on the chromosome of interest (if provided)
    let mut a2 = if args.chrom != "all" {
        a.into_iter()
            .filter(|r| r.reference_sequence_name() == args.chrom)
            .collect()
    } else {
        a
    };

    //  Pass2. Intersect out (subtract non-desired regions) (if provided)
    let a3 = match &args.bed_out {
        Some(bed_out) => {
            let b = load_bed_records(&bed_out)?;
            a2.retain(|r| !b.contains(r));
            a2
        }
        None => a2,
    };

    // Pass3. Makewindows in final bed file based on the window sizes specified as argument
    let final_bed = generate_windows_for_records(a3, args.bin)?;
    Ok(final_bed)
}

fn generate_windows_for_records(
    record: Vec<Record<3>>,
    window_size: usize,
) -> Result<Vec<Record<3>>> {
    let mut all_records = vec![];
    for r in record {
        let records = generate_windows_for_record(r, window_size)?;
        all_records.extend(records);
    }

    Ok(all_records)
}

fn generate_windows_for_record(record: Record<3>, window_size: usize) -> Result<Vec<Record<3>>> {
    let mut windows = Vec::new();

    let start = record.start_position().get();
    let end = record.end_position().get();

    for window_start in (start..end).step_by(window_size) {
        let window_end = cmp::min(window_start + window_size, end);
        let window = bed::Record::<3>::builder()
            .set_reference_sequence_name(record.reference_sequence_name())
            .set_start_position(window_start.try_into()?)
            .set_end_position(window_end.try_into()?)
            .build()?;

        windows.push(window);
    }

    Ok(windows)
}

fn get_records(args: &Args) -> Result<Vec<bed::Record<3>>> {
    if let Some(bed) = &args.bed {
        load_bed_records(&bed)
    } else {
        // We get the bed file based on all coordenates from reference fasta file
        let in_fasta = faidx::Reader::from_path("./test/chr10.fa")?;

        let mut list: Vec<bed::Record<3>> = vec![];
        for i in 0..in_fasta.n_seqs() {
            let name = in_fasta.seq_name(i as i32)?; // TODO review cast

            let start_position = 1;
            let end_position = in_fasta.fetch_seq_len(&name);

            let record = bed::Record::<3>::builder()
                .set_reference_sequence_name(name)
                .set_start_position(start_position.try_into()?)
                .set_end_position((end_position as usize).try_into()?)
                .build()?;

            list.push(record);
        }
        Ok(list)
    }
}

fn load_bed_records(path: &str) -> Result<Vec<bed::Record<3>>> {
    let mut reader = File::open(path).map(BufReader::new).map(bed::Reader::new)?;

    let records_iter = reader.records::<3>();
    let mut records: Vec<bed::Record<3>> = vec![];
    for record_res in records_iter {
        records.push(record_res?);
    }
    Ok(records)
}

#[cfg(test)]
mod test {
    use super::{make_windows, Args};

    #[test]
    fn test_make_windows() {
        // corresponds to example at https://github.com/cortes-ciriano-lab/SComatic/blob/main/docs/SComaticExample.md#step-2-collecting-base-count-information
        let args = Args {
            fasta: "./test/Example.scrnaseq.bam".to_string(),
            chrom: "all".to_string(),
            bed: None,
            bed_out: None,
            bin: 50000,
        };

        let records = make_windows(&args).unwrap();

        // println!("records: {:?}", records);
        assert_eq!(records.len(), 2676);

        // for now test just the 2 first and 2 last elements
        let first = &records[0];
        assert_eq!(first.reference_sequence_name(), "chr10");
        assert_eq!(first.start_position().get(), 1);
        assert_eq!(first.end_position().get(), 50001);
        let second = &records[1];
        assert_eq!(second.reference_sequence_name(), "chr10");
        assert_eq!(second.start_position().get(), 50001);
        assert_eq!(second.end_position().get(), 100001);
        let second_to_last = &records[records.len() - 2];
        assert_eq!(second_to_last.reference_sequence_name(), "chr10");
        assert_eq!(second_to_last.start_position().get(), 133700001);
        assert_eq!(second_to_last.end_position().get(), 133750001);
        let last = &records[records.len() - 1];
        assert_eq!(last.reference_sequence_name(), "chr10");
        assert_eq!(last.start_position().get(), 133750001);
        assert_eq!(last.end_position().get(), 133797422);
    }
}
