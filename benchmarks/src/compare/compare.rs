use glob::glob;
use needletail::{Sequence, parse_fastx_file};
use sassy::profiles::Iupac;
use sassy::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::Command;
use std::time::Instant;

const MAX_FLANK_EDITS: usize = 15;
const FLANK_SEQ: &[u8] = b"GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA";

pub struct Dorado {
    pub exec_path: String,
}

pub struct Barbell {
    pub exec_path: String,
}

pub struct Flexiplex {
    pub exec_path: String,
}

trait Tool {
    fn new(exec_path: &str) -> Self;
    // Other params we considered fixed and edit in this source file
    fn run(
        &self,
        fastq_file: &str,
        output_folder: &str,
        threads: usize,
        extra_file: Option<String>,
    ) -> Result<(), Box<dyn std::error::Error>>;
    /// Parse output expecting a file like read_id > barcode_id, tsv
    fn parse_output(
        &self,
        output_folder: &str,
        anno_out_file: &str,
        trimmed_out_file: &str,
        extra_file: Option<String>,
    );
}

impl Tool for Dorado {
    fn new(exec_path: &str) -> Self {
        // check if exec path exists
        if !Path::new(&exec_path).exists() {
            panic!("Dorado executable not found at {exec_path}");
        }
        Self {
            exec_path: exec_path.to_string(),
        }
    }

    fn run(
        &self,
        fastq_file: &str,
        output_folder: &str,
        threads: usize,
        _extra_file: Option<String>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        // Create output folder if not exists
        if !Path::new(output_folder).exists() {
            std::fs::create_dir_all(output_folder).expect("Failed to create output folder");
        }

        let cmd = format!(
            "{0} demux --kit-name SQK-RBK114-96 -o {output_folder} --emit-fastq {fastq_file} -t {threads} --min-score 0.2 --min-score-diff 0.1",
            self.exec_path
        );
        let result = Command::new("bash").arg("-c").arg(cmd).output()?;

        if !result.status.success() {
            return Err(format!(
                "Command failed: {}",
                String::from_utf8_lossy(&result.stderr)
            )
            .into());
        }
        Ok(())
    }

    fn parse_output(
        &self,
        output_folder: &str,
        anno_out_file: &str,
        trimmed_out_file: &str,
        extra_file: Option<String>,
    ) {
        // We have to go over each fastq file in the output folder,
        // and read it using needletail, extracting barcode id from file name
        // and then storing read_id > barcode_id
        let mut searcher = Searcher::<Iupac>::new_rc();

        let output_file_handle = File::create(anno_out_file).expect("Failed to create output file");
        let mut writer = BufWriter::new(output_file_handle);

        let trimmed_out_handle =
            File::create(trimmed_out_file).expect("Failed to create trim output file");
        let mut trimmed_writer = BufWriter::new(trimmed_out_handle);

        for entry in glob(&format!("{output_folder}/*.fastq")).expect("Failed to read glob pattern")
        {
            let fastq_path = entry.unwrap();

            // Split path on / and get last element
            let barcode_id = fastq_path
                .file_name()
                .unwrap()
                .to_str()
                .unwrap()
                .strip_suffix(".fastq")
                .unwrap()
                .split("_")
                .last()
                .unwrap();

            // We skip unclassified reads
            if barcode_id == "unclassified" {
                continue;
            }

            let mut reader = parse_fastx_file(&fastq_path).expect("valid path/file");
            while let Some(record) = reader.next() {
                let record = record.unwrap();
                let norm_seq = record.seq().normalize(false).to_vec();
                let seq = String::from_utf8_lossy(&norm_seq);
                let read_id = String::from_utf8_lossy(record.id());
                let seq_len = record.seq().len();
                let flank_matches = searcher.search(FLANK_SEQ, &norm_seq, MAX_FLANK_EDITS);
                let n_flank_matches = flank_matches.len();
                let anno_line = format!("{read_id}\t{barcode_id}\t{seq_len}\t{n_flank_matches}\n");
                let seq_line = format!(">{read_id}\n{seq}\n");

                writer
                    .write_all(anno_line.as_bytes())
                    .expect("Failed to write annotation");

                trimmed_writer
                    .write_all(seq_line.as_bytes())
                    .expect("Failed to write annotation");
            }
        }
    }
}

impl Tool for Barbell {
    fn new(exec_path: &str) -> Self {
        // check if exec path exists
        if !Path::new(&exec_path).exists() {
            panic!("Barbell executable not found at {exec_path}");
        }
        Self {
            exec_path: exec_path.to_string(),
        }
    }

    fn run(
        &self,
        fastq_file: &str,
        output_folder: &str,
        threads: usize,
        extra_file: Option<String>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        // Create output dir if not existing
        if !Path::new(output_folder).exists() {
            std::fs::create_dir_all(output_folder).expect("Failed to create output folder");
        }

        let failed_out = format!("{output_folder}/failed_trimmed.txt");

        let cmd = format!(
            "{0} kit -k \"SQK-RBK114-96\" -i {fastq_file} -o {output_folder} -t {threads} --failed-out {failed_out} --maximize",
            self.exec_path
        );
        let result: std::process::Output = Command::new("bash").arg("-c").arg(cmd).output()?;

        if !result.status.success() {
            return Err(format!(
                "Command failed: {}",
                String::from_utf8_lossy(&result.stderr)
            )
            .into());
        }
        Ok(())
    }

    fn parse_output(
        &self,
        output_folder: &str,
        anno_out_file: &str,
        trimmed_out_file: &str,
        extra_file: Option<String>,
    ) {
        // We have to go over each fastq file in the output folder,
        // and read it using needletail, extracting barcode id from file name
        // and then storing read_id > barcode_id
        let mut searcher = Searcher::<Iupac>::new_rc();

        let output_file_handle = File::create(anno_out_file).expect("Failed to create output file");
        let mut writer = BufWriter::new(output_file_handle);

        let trimmed_out_handle =
            File::create(trimmed_out_file).expect("Failed to create trim output file");
        let mut trimmed_writer = BufWriter::new(trimmed_out_handle);

        for entry in glob(&format!("{output_folder}/*.fastq")).expect("Failed to read glob pattern")
        {
            let fastq_path = entry.unwrap();
            //println!("Fastq file: {}", fastq_path.display());
            // Split path on / and get last element
            let barcode_id = fastq_path
                .file_name()
                .unwrap()
                .to_str()
                .unwrap()
                .split("_")
                .next()
                .unwrap();

            let mut reader = parse_fastx_file(&fastq_path).expect("valid path/file");
            while let Some(record) = reader.next() {
                let record = record.unwrap();
                let norm_seq = record.seq().normalize(false).to_vec();
                let seq = String::from_utf8_lossy(&norm_seq);
                let read_id = String::from_utf8_lossy(record.id());
                let flank_matches = searcher.search(FLANK_SEQ, &norm_seq, MAX_FLANK_EDITS);
                let n_flank_matches = flank_matches.len();
                let seq_len = record.seq().len();
                let anno_line = format!("{read_id}\t{barcode_id}\t{seq_len}\t{n_flank_matches}\n");
                let seq_line = format!(">{read_id}\n{seq}\n");

                writer
                    .write_all(anno_line.as_bytes())
                    .expect("Failed to write annotation");

                trimmed_writer
                    .write_all(seq_line.as_bytes())
                    .expect("Failed to write annotation");
            }
        }
    }
}

impl Tool for Flexiplex {
    fn new(exec_path: &str) -> Self {
        // check if exec path exists
        if !Path::new(&exec_path).exists() {
            panic!("Flexiplex executable not found at {exec_path}");
        }
        Self {
            exec_path: exec_path.to_string(),
        }
    }

    fn run(
        &self,
        fastq_file: &str,
        output_folder: &str,
        threads: usize,
        extra_file: Option<String>,
    ) -> Result<(), Box<dyn std::error::Error>> {
        // Create output dir if not existing
        if !Path::new(output_folder).exists() {
            std::fs::create_dir_all(output_folder).expect("Failed to create output folder");
        }
        let barcode_file = extra_file.expect("Flexiplex needs a barcode file");

        let output_file = format!("{output_folder}/classified_reads.fastq");
        let output_handle = File::create(&output_file)?;

        let cmd = format!(
            "{0} -x GCTTGGGTGTTTAACC -b ???????????????????????? -x GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA -e 7 -f 20 -p {threads} -k {barcode_file} -s true -n {output_folder} {fastq_file}",
            self.exec_path
        );

        let mut child = Command::new("bash")
            .arg("-c")
            .arg(cmd)
            .stdout(output_handle)
            .spawn()?;

        let status = child.wait()?;
        if !status.success() {
            return Err("Flexiplex command failed".into());
        }
        Ok(())
    }

    fn parse_output(
        &self,
        output_folder: &str,
        anno_out_file: &str,
        trimmed_out_file: &str,
        extra_file: Option<String>,
    ) {
        let mut searcher = Searcher::<Iupac>::new_rc();

        // Flexiplex outputs to a single file: {output_folder}/classified_reads.fastq
        let flexiplex_output = format!("{output_folder}/classified_reads.fastq");

        if !Path::new(&flexiplex_output).exists() {
            eprintln!("Flexiplex output file not found: {}", flexiplex_output);
            return;
        }

        // Create mapping of sequence to barcode id using extra file
        if extra_file.is_none() {
            eprint!("We need barcode file to link sequence to barcode name");
            return;
        }

        // We read line by line, creat mapping seq > barcode
        let handle = File::open(extra_file.unwrap()).expect("barcode file not found");
        let reader = BufReader::new(handle);
        let mut barcode_map = HashMap::new();
        for line in reader.lines() {
            let line = line.expect("Failed to read line from barcode file");
            let parts: Vec<String> = line.split('\t').map(|s| s.to_owned()).collect();
            if parts.len() >= 2 {
                barcode_map.insert(parts[0].clone(), parts[1].clone());
            }
        }

        let output_file_handle = File::create(anno_out_file).expect("Failed to create output file");
        let mut writer = BufWriter::new(output_file_handle);

        let trimmed_out_handle =
            File::create(trimmed_out_file).expect("Failed to create trim output file");
        let mut trimmed_writer = BufWriter::new(trimmed_out_handle);

        println!("output file: {}", flexiplex_output);
        let mut reader = parse_fastx_file(&flexiplex_output).expect("valid path/file");
        while let Some(record) = reader.next() {
            let record = record.unwrap();
            let norm_seq = record.seq().normalize(false).to_vec();
            let seq = String::from_utf8_lossy(&norm_seq);
            let read_id = String::from_utf8_lossy(record.id());
            // println!("Read id: {}", read_id);

            // Parse it
            let barcode_seq = read_id
                // .strip_prefix("@")
                // .unwrap()
                .split("_")
                .next()
                .unwrap();
            let barcode_id = barcode_map.get(barcode_seq).unwrap();
            // CTTTCGTTGTTGACTCGACGGTAG_#2b86d04f-daf8-48a2-a083-b0b50b1138d6_-1of1
            let read_id = read_id
                .split("#")
                .nth(1)
                .unwrap()
                .split("_")
                .next()
                .unwrap();
            // TODO: we could use the 1 of 1 for double split cases to make the ids unique
            // but will make joining later more complex
            let flank_matches = searcher.search(FLANK_SEQ, &norm_seq, MAX_FLANK_EDITS);
            let n_flank_matches = flank_matches.len();
            let seq_len = record.seq().len();
            let anno_line = format!("{read_id}\t{barcode_id}\t{seq_len}\t{n_flank_matches}\n");
            let seq_line = format!(">{read_id}\n{seq}\n");

            writer
                .write_all(anno_line.as_bytes())
                .expect("Failed to write annotation");

            trimmed_writer
                .write_all(seq_line.as_bytes())
                .expect("Failed to write trimmed sequence");
        }
    }
}

pub fn run_all_tools(
    fastq_file: &str,
    output_folder: &str,
    threads: usize,
    extra_file: Option<String>,
    dorado_exec_path: &str,
    barbell_exec_path: &str,
    flexiplex_exec_path: &str,
) {
    // We create additional output folders for each of the tools
    let output_folder = format!("{output_folder}/all_tools");
    let dorado_output_folder = format!("{output_folder}/dorado");
    let barbell_output_folder = format!("{output_folder}/barbell");
    let flexiplex_output_folder = format!("{output_folder}/flexiplex");
    let annotation_output_folder = format!("{output_folder}/annotation");
    let trimmed_output_folder = format!("{output_folder}/trimmed");

    // Create annotationa nd trimmed folders if not existing
    std::fs::create_dir_all(&annotation_output_folder)
        .expect("Failed to create annotation output folder");
    std::fs::create_dir_all(&trimmed_output_folder)
        .expect("Failed to create trimmed output folder");

    // We first create a file that has the original (untrimmed) sequence length
    // for each of the reads
    let untrimmed_lengths_file = format!("{output_folder}/untrimmed_lengths.tsv");
    let untrimmed_lengths_handle =
        File::create(untrimmed_lengths_file).expect("Failed to create untrimmed lengths file");
    let mut untrimmed_lengths_writer = BufWriter::new(untrimmed_lengths_handle);
    let mut reader = parse_fastx_file(fastq_file).expect("valid path/file");
    while let Some(record) = reader.next() {
        let record = record.unwrap();
        let seq_len = record.seq().len();
        let read_id = String::from_utf8_lossy(record.id());
        let read_len_line = format!("{read_id}\t{seq_len}\n");
        untrimmed_lengths_writer
            .write_all(read_len_line.as_bytes())
            .expect("Failed to write untrimmed lengths");
    }

    // // -- Dorado --
    // let dorado = Dorado::new(dorado_exec_path);
    // let start_time = Instant::now();
    // println!("Running Dorado");
    // dorado
    //     .run(fastq_file, &dorado_output_folder, threads, None)
    //     .unwrap();
    // let dorado_time: std::time::Duration = start_time.elapsed();
    // println!("Dorado time: {:?}", dorado_time);
    // dorado.parse_output(
    //     &dorado_output_folder,
    //     &format!("{annotation_output_folder}/dorado_parsed.tsv"),
    //     &format!("{trimmed_output_folder}/dorado_trimmed.fasta"),
    //     None,
    // );

    // -- Barbell --
    println!("Running Barbell");
    let barbell = Barbell::new(barbell_exec_path);
    let start_time = Instant::now();
    barbell
        .run(fastq_file, &barbell_output_folder, threads, None)
        .unwrap();
    let barbell_time = start_time.elapsed();
    println!("Barbell time: {:?}", barbell_time);
    barbell.parse_output(
        &barbell_output_folder,
        &format!("{annotation_output_folder}/barbell_parsed.tsv"),
        &format!("{trimmed_output_folder}/barbell_trimmed.fasta"),
        None,
    );

    // // // -- Flexiplex --
    // let flexiplex = Flexiplex::new(flexiplex_exec_path);
    // let start_time = Instant::now();
    // println!("Running Flexiplex");
    // flexiplex
    //     .run(
    //         fastq_file,
    //         &flexiplex_output_folder,
    //         threads,
    //         extra_file.clone(),
    //     )
    //     .unwrap();
    // let flexiplex_time = start_time.elapsed();
    // println!("Flexiplex time: {:?}", flexiplex_time);
    // flexiplex.parse_output(
    //     &flexiplex_output_folder,
    //     &format!("{annotation_output_folder}/flexiplex_parsed.tsv"),
    //     &format!("{trimmed_output_folder}/flexiplex_trimmed.fasta"),
    //     extra_file.clone(),
    // );
    // println!("All done!");
    // println!("Timings");
    // println!("Dorado: {:?}", dorado_time);
    // println!("Barbell: {:?}", barbell_time);
    // println!("Flexiplex: {:?}", flexiplex_time);
}

/*
-x GCTTGGGTGTTTAACC -b ???????????????????????? "
        "-x GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA -e 6 -f 17 -p {threads} -k {input[1]} "
        "-s true -n {wildcards.group} {input[0]}
         */

#[cfg(test)]
mod tests {
    use super::*;
}
