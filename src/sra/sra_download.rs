use anyhow::Context;
use std::io::Write;
use tempfile::NamedTempFile;

/// Download an SRA accession and write FASTQ data to a temp file.
pub fn download_sra_to_tempfile(accession: &str, threads: usize) -> anyhow::Result<NamedTempFile> {
    let mut opts = xsra::DumpOptions::new(accession);
    opts.threads = threads as u64;
    let bytes = xsra::dump_with_options(opts)
        .with_context(|| format!("Failed to download SRA accession '{accession}' with xsra"))?;
    let mut file = NamedTempFile::with_suffix(".fastq")
        .context("Failed to create temporary FASTQ file for SRA download")?;
    file.write_all(&bytes)
        .context("Failed to write downloaded SRA FASTQ to temp file")?;
    file.flush()
        .context("Failed to flush downloaded SRA FASTQ temp file")?;
    Ok(file)
}
