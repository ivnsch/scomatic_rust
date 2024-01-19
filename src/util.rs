use anyhow::Result;
use rust_htslib::{bam, bam::Read};

/// debug
pub fn to_sam(bam_path: &str, sam_path: &str) -> Result<()> {
    let mut bam_reader = bam::Reader::from_path(bam_path)?;

    // Create a SAM header based on the BAM header
    let header = bam_reader.header().to_owned();
    let sam_header = bam::Header::from_template(&header);

    let mut sam_writer = bam::Writer::from_path(sam_path, &sam_header, bam::Format::Sam)
        .expect("Error creating SAM writer");

    for record in bam_reader.records() {
        let record = record.expect("Error reading BAM record");
        sam_writer.write(&record).expect("Error writing SAM record");
    }

    Ok(())
}
