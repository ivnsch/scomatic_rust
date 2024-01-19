use anyhow::{anyhow, Result};
use clap::{arg, command, Parser};
use polars::prelude::*;
use polars::series::SeriesIter;
use polars_lazy::dsl::GetOutput;
use rust_htslib::bam::index::Type;
use rust_htslib::bam::record::{Aux, Cigar};
use rust_htslib::bam::{Header, Record};
use rust_htslib::{bam, bam::Read};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    #[arg(long, help = "BAM file to be analysed (Sorted by coordinate)")]
    bam: String,
    #[arg(
        long = "meta",
        help = "Metadata file mapping cell barcodes to cell type information"
    )]
    txt: String,
    #[arg(long, help = "Out directory")]
    outdir: String,
    #[arg(long = "id", help = "Sample ID")]
    donor: String,
    #[arg(long, help = "BAM file to be analysed (Sorted by coordinate)")]
    tissue: Option<String>,
    #[arg(
        long = "max_nM",
        help = "Maximum number of mismatches permitted to consider reads for analysis. By default, this filter is switched off, although we recommed using --max_nM 5. If applied, this filter requires having the nM tag in the bam file. [Default: Switched off]"
    )]
    max_nm: Option<i8>,
    #[arg(
        long = "max_NH",
        help = "Maximum number of alignment hits permitted to consider reads for analysis. By default, this filter is switched off, although we recommend using --max_NH 1. This filter requires having the NH tag in the bam file. [Default: Switched off]"
    )]
    max_nh: Option<i8>,
    #[arg(
        long = "min_MQ",
        help = "Minimum mapping quality required to consider reads for analysis. Set this value to 0 to switch this filter off. --min_MQ 255 is recommended for RNA data, and --min_MQ 30 for DNA data. [Default: 255]",
        default_value_t = 255
    )]
    min_mapq: u32,
    #[arg(
        long = "n_trim",
        help = "Number of bases trimmed by setting the base quality to 0 at the beginning and end of each read [Default: 0]",
        default_value_t = 0
    )]
    n_trim: u32,
}

pub fn main() -> Result<()> {
    let args = &Args::parse();

    let start = Instant::now();

    let annotations = load_annotations(&args.txt, args.tissue.clone())?;

    let mut stats = split_bam(args, &annotations)?;

    let elapsed = start.elapsed();
    stats.insert("Total_time".to_string(), elapsed.as_millis());

    write_report(args, stats)?;

    index_bam_files(args, &annotations.cell_types)?;

    // debugging
    // to_sam(args.bam, "./sam_export.sam");

    Ok(())
}

fn split_bam(args: &Args, annotations: &Annotations) -> Result<HashMap<String, u128>> {
    if annotations.is_empty() {
        return Err(anyhow!(
            "Warning: No cell barcodes found in the --meta file"
        ));
    }

    let mut infile = bam::Reader::from_path(&args.bam)?;

    // Create and open out bam files
    let header = bam::Header::from_template(infile.header());
    let mut dict_files = create_cell_type_files(args, &annotations.cell_types, header)?;

    // Start read counts
    let mut total_reads = 0;
    let mut stats = HashMap::from([
        ("Total_reads".to_string(), 0),
        ("Pass_reads".to_string(), 0),
        ("CB_not_found".to_string(), 0),
        ("CB_not_matched".to_string(), 0),
    ]);

    // Check reads and split them in bam files

    for record in infile.records() {
        let mut record: Record = record?;

        // println!("record: {:?}", String::from_utf8_lossy(read.qname()));

        let cell_type = get_cell_type_and_update_stats(
            args,
            &record,
            &mut stats,
            &mut total_reads,
            annotations,
        );

        // Overwrite qualities
        if args.n_trim > 0 {
            let qualities: Vec<u8> = get_qualities(args.n_trim, &record)?;
            set_qualities(&mut record, &qualities);
        }

        // Write read
        if let Some(cell_type) = cell_type {
            if let Some(writer) = dict_files.get_mut(&cell_type) {
                writer.write(&record)?;
            } else {
                println!("Key not found: {}", cell_type);
            }
        }
    }

    Ok(stats)
}

fn get_cell_type_and_update_stats(
    args: &Args,
    record: &Record,
    stats: &mut HashMap<String, u128>,
    total_reads: &mut i32,
    annotations: &Annotations,
) -> Option<String> {
    stats
        .entry("Total_reads".to_string())
        .and_modify(|c| *c += 1);
    *total_reads += 1;

    if *total_reads % 5000000 == 0 {
        println!("Number of records already processed: {}", total_reads);
    }

    if let Some(cell_type) = get_cell_type_with_barcode_and_update_stats(record, stats, annotations)
    {
        // Final filters
        match get_filter_str(args, record) {
            Some(s) => {
                *stats.entry(s).or_insert(0) += 1;
            }
            None => {
                stats
                    .entry("Pass_reads".to_string())
                    .and_modify(|c| *c += 1);
            }
        }
        Some(cell_type)
    } else {
        None
    }
}

fn get_cell_type_with_barcode_and_update_stats(
    record: &Record,
    stats: &mut HashMap<String, u128>,
    annotations: &Annotations,
) -> Option<String> {
    match get_barcode(record) {
        Ok(barcode) => match annotations.index.get(&barcode) {
            Some(c) => return Some(c.to_string()),
            None => {
                stats
                    .entry("CB_not_matched".to_string())
                    .and_modify(|c| *c += 1);
                return None;
            }
        },
        Err(_) => {
            stats
                .entry("CB_not_found".to_string())
                .and_modify(|c| *c += 1);
            return None;
        }
    }
}

fn load_annotations(txt: &str, tissue: Option<String>) -> Result<Annotations> {
    let df = CsvReader::from_path(txt)
        .unwrap()
        .with_separator(*"\t".as_bytes().first().unwrap())
        .has_header(true)
        .finish()?;

    // Clean index column
    let mut out = df
        .lazy()
        .with_columns([col("Index")
            .str()
            .replace_all(lit("[-|.|*|$]"), lit(""), false)
            .alias("Index_clean")])
        .with_columns([col("Cell_type")
            .str()
            .replace_all(lit(" "), lit(""), false)
            .alias("Cell_type_clean")]);

    // If tissue provided, append tissue id to cell type to be recognised in downstream analysis
    if let Some(tissue) = tissue {
        out = append_tissue_to_cell_type(&tissue, out)?;
    }

    dataframe_to_annotations(out.collect()?)
}

fn append_tissue_to_cell_type(tissue: &str, out: LazyFrame) -> Result<LazyFrame> {
    let tissue = tissue.replace(' ', "");

    fn map_value(cell_type: Option<&str>, tissue: &str) -> String {
        match cell_type {
            Some(s) => format!("{}__{}", tissue, s),
            None => "".into(),
        }
    }

    fn map_series(e: Series, tissue: &str) -> Result<Option<Series>, PolarsError> {
        Ok(Some(
            e.str()?
                .into_iter()
                .map(|cell_type| map_value(cell_type, tissue))
                .collect(),
        ))
    }

    Ok(out.with_columns([col("Cell_type_clean").map(
        move |series| map_series(series, &tissue.clone()),
        GetOutput::from_type(DataType::String),
    )]))
}

struct Annotations {
    /// barcode -> cell type
    index: HashMap<String, String>,
    /// cell types from index, for convenient access
    cell_types: HashSet<String>,
}

impl Annotations {
    fn is_empty(&self) -> bool {
        // could also be cell_types
        self.index.is_empty()
    }
}

// https://stackoverflow.com/questions/72440403/iterate-over-rows-polars-rust/72443329#72443329
fn dataframe_to_annotations(df: DataFrame) -> Result<Annotations> {
    let mut dict: HashMap<String, String> = HashMap::new();
    let mut cell_types = HashSet::new();

    let mut iters = df
        .columns(["Index_clean", "Cell_type_clean"])?
        .iter()
        .map(|s| s.iter())
        .collect::<Vec<_>>();

    fn next_str(series_iter: &mut SeriesIter<'_>) -> String {
        series_iter
            .next()
            .expect("should have as many iterations as rows")
            .get_str()
            .expect("should be a string")
            .to_string()
    }

    for _row in 0..df.height() {
        let iter = &mut iters.iter_mut();

        let key_iter = iter.next().expect("should be set");
        let value_iter = iter.next().expect("should be set");

        let index = next_str(key_iter);
        let cell_type = next_str(value_iter);

        dict.insert(index, cell_type.clone());
        cell_types.insert(cell_type);
    }

    Ok(Annotations {
        index: dict,
        cell_types,
    })
}

/// Creates a file for each cell type and returns a hashmap (cell type -> file writer)
fn create_cell_type_files(
    args: &Args,
    all_cell_types: &HashSet<String>,
    header: Header,
) -> Result<HashMap<String, bam::Writer>> {
    let mut dict_files = HashMap::<String, bam::Writer>::new();

    for cell_type in all_cell_types {
        let path = format!("{}/{}.{}.bam", args.outdir, args.donor, cell_type);
        File::create(path.clone())?;
        let writer = bam::Writer::from_path(path, &header, bam::Format::Bam)?;
        dict_files.insert(cell_type.to_string(), writer);
    }

    Ok(dict_files)
}

fn set_qualities(record: &mut Record, qualities: &[u8]) {
    let record_clone = record.clone();
    let qname = record_clone.qname();
    let seq = record_clone.seq();
    let cigar = record_clone.cigar().take();
    // Substitute qualities
    record.set(qname, Some(&cigar), &seq.as_bytes(), qualities);
}

/// Only for PASS reads
/// Trim last and first bases of the read (reduce quality) if specified
/// It does not consider the soft-clip bases at the beginning/end of the read for the counting
/// It also considers the expected adapter lengths (up to 30) of 10x library prep to remove bases in long soft-clip sequences
fn get_qualities(n_trim: u32, record: &Record) -> Result<Vec<u8>> {
    let cigars = record.cigar();

    let mut trim_start = n_trim;
    let mut trim_end = n_trim;

    if cigars.len() > 1 {
        // Check the number of bases to trim at the beginning of the record
        // unwrap: we checked for > 1
        if let Some(trim) = trim_if_softclip(cigars.first().unwrap(), n_trim) {
            trim_start = trim
        }
        // Check the number of bases to trim at the end of the record
        // unwrap: we checked for > 1
        if let Some(trim) = trim_if_softclip(cigars.last().unwrap(), n_trim) {
            trim_end = trim
        }
    }

    // Set base qualities to 0 at the beginning and the end of the record if specified
    let mut qualities = record.qual().to_vec();
    let len = qualities.len();
    let start_range = 0..(trim_start as usize);
    let end_range = (len - trim_end as usize)..len;
    let indexes = start_range.chain(end_range);
    for i in indexes {
        qualities[i] = 0;
    }

    Ok(qualities)
}

fn trim_if_softclip(cigar: &Cigar, n_trim: u32) -> Option<u32> {
    match cigar {
        // Extension of the soft-clipped bases
        bam::record::Cigar::SoftClip(n) => {
            // As there are some library preparation adapters inside the 10x prep protocol that
            // are not properly removed (up to 30 bps), if we observed that the softclip bases are around more than 20, we
            // prefer to be conservative and remove the expected 30 pbs + the n_trim parameter
            Some(if (&20..&30).contains(&n) {
                30 + n_trim
            } else {
                n + n_trim
            })
        }
        _ => None,
    }
}

#[cfg(test)]
mod test {
    use super::{load_annotations, split_bam, Args};

    #[test]
    fn test_generates_expected_report() {
        // corresponds to example at https://github.com/cortes-ciriano-lab/SComatic/blob/main/docs/SComaticExample.md#step-1-splitting-alignment-file-in-cell-type-specific-bams
        // note 255 is default, but need to pass explicitly with struct
        let args = Args {
            bam: "./test/Example.scrnaseq.bam".to_string(),
            txt: "./test/Example.cell_barcode_annotations.tsv".to_string(),
            outdir: ".".to_string(),
            donor: "Example".to_string(),
            tissue: None,
            max_nm: Some(5),
            max_nh: Some(1),
            min_mapq: 255,
            n_trim: 5,
        };

        let annotations = load_annotations(&args.txt, args.tissue.clone()).unwrap();

        let stats_res = split_bam(&args, &annotations);
        assert!(stats_res.is_ok());

        let stats = stats_res.unwrap();

        assert_eq!(*stats.get("nM").unwrap(), 32);
        assert_eq!(*stats.get("NH;MAPQ").unwrap(), 11);
        assert_eq!(*stats.get("Pass_reads").unwrap(), 459);
        assert_eq!(*stats.get("CB_not_found").unwrap(), 10);
        assert_eq!(*stats.get("Total_reads").unwrap(), 1158);
        assert_eq!(*stats.get("nM;NH;MAPQ").unwrap(), 5);
        assert_eq!(*stats.get("CB_not_matched").unwrap(), 641);
    }
}

fn get_barcode(record: &Record) -> Result<String> {
    match get_aux_str(record, b"CB") {
        Ok(s) => {
            let parts: Vec<&str> = s.split('-').collect::<Vec<&str>>();
            match parts.first() {
                Some(part) => Ok(part.to_string()),
                None => Err(anyhow!("Unexpected format: {:?}", parts)),
            }
        }
        Err(e) => Err(e),
    }
}

fn get_aux_str(record: &Record, tag: &[u8]) -> Result<String> {
    match record.aux(tag) {
        Ok(aux) => match aux {
            Aux::String(s) => Ok(s.to_string()),
            _ => panic!("Aux was not a string: {:?}", aux),
        },
        Err(e) => Err(e.into()),
    }
}

fn get_i8_str(record: &Record, tag: &[u8]) -> Result<i8> {
    match record.aux(tag) {
        Ok(aux) => match aux {
            Aux::I8(i) => Ok(i),
            _ => panic!("Aux was not an i8: {:?}", aux),
        },
        Err(e) => Err(e.into()),
    }
}

fn get_filter(args: &Args, record: &Record) -> Vec<String> {
    let mut filter: Vec<String> = vec![];

    // Check number of mismatches in the record
    if let Some(max_nm) = args.max_nm {
        match get_i8_str(record, b"nM") {
            Ok(i) => {
                if i > max_nm {
                    filter.push("nM".to_string())
                }
            }
            Err(_) => filter.push("nM_not_found".to_owned()),
        }
    }

    // Check number of hits
    if let Some(max_nh) = args.max_nh {
        match record.aux(b"NH") {
            Ok(aux) => match aux {
                Aux::I8(a) => {
                    if a > max_nh {
                        filter.push("NH".to_string())
                    }
                }
                _ => {
                    panic!("unexpected data type: {:?}", aux)
                }
            },
            Err(_) => {
                filter.push("NH_not_found".to_owned());
            }
        }
    }

    // Check mapping quality
    if args.min_mapq > 0 && record.mapq() < args.min_mapq as u8 {
        filter.push("MAPQ".to_owned());
        // except:
        // 		FILTER.append('MAPQ_not_found')
        // no exception here - library always returns a vamlue for mapq
    }

    filter
}

// Final filters
fn get_filter_str(args: &Args, record: &Record) -> Option<String> {
    let filter: Vec<String> = get_filter(args, record);

    if !filter.is_empty() {
        Some(filter.join(";").to_string())
    } else {
        None
    }
}

fn write_report(args: &Args, stats: HashMap<String, u128>) -> Result<()> {
    let outfile_report: String = format!("{}/{}.report.txt", args.outdir, args.donor);

    let (keys, values): (Vec<String>, Vec<String>) = stats
        .iter()
        .map(|(k, v)| (k.to_string().clone(), v.to_string().clone()))
        .unzip();

    let keys_series = Series::new("keys", &keys);
    let values_series = Series::new("values", &values);

    let mut data_df = DataFrame::new(vec![keys_series, values_series])?;

    let mut file = File::create(outfile_report)?;
    CsvWriter::new(&mut file).finish(&mut data_df)?;

    Ok(())
}

fn index_bam_files(args: &Args, all_cell_types: &HashSet<String>) -> Result<()> {
    for cell_type in all_cell_types {
        let path = format!("{}/{}.{}.bam", args.outdir, args.donor, cell_type);
        File::open(path.clone()).expect("Failed to create file");
        bam::index::build(&path, None, Type::Bai, 1)?
    }

    Ok(())
}
