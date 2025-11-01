use barbell::annotate::barcodes::{BarcodeGroup, BarcodeType};
use barbell::annotate::edit_model::get_edit_cut_off;
use barbell::annotate::search_strategies::{SearchMode, SearchStrategy};
use barbell::annotate::searcher::Demuxer;
use needletail::parse_fastx_file;
use std::env;
use std::path::PathBuf;

fn main() -> anyhow::Result<()> {
    // Get command line args for flexibility
    let args: Vec<String> = env::args().collect();

    // Default paths - can be overridden
    let fastq_path = if args.len() > 1 {
        args[1].clone()
    } else {
        // Try to find a test fastq file
        let candidates = vec!["tests/sample.fastq", "test.fastq", "combined.fastq"];

        candidates
            .into_iter()
            .find(|p| PathBuf::from(p).exists())
            .expect("No FASTQ file found. Please provide path as first argument.")
            .to_string()
    };

    let kit = if args.len() > 2 {
        args[2].clone()
    } else {
        "SQK-RBK114-96".to_string()
    };

    let iterations = if args.len() > 3 {
        args[3].parse::<usize>().unwrap_or(10)
    } else {
        10 // Run 10 times through the data by default
    };

    println!("Profiling demux with:");
    println!("  FASTQ: {}", fastq_path);
    println!("  Kit: {}", kit);
    println!("  Iterations: {}", iterations);
    println!();

    // Setup query groups from kit
    let query_groups: Vec<BarcodeGroup> = BarcodeGroup::new_from_kit(&kit, false);
    let query_groups: Vec<BarcodeGroup> = query_groups
        .into_iter()
        .map(|mut query_group| {
            let edit_cut_off = get_edit_cut_off(query_group.get_effective_len());
            if query_group.barcode_type == BarcodeType::Ftag
                || query_group.barcode_type == BarcodeType::Rtag
            {
                let barcode_cutoff = get_edit_cut_off(query_group.get_barcode_len());
                query_group.set_barcode_threshold(barcode_cutoff);
            }
            query_group.set_flank_threshold(edit_cut_off);
            query_group
        })
        .collect();

    // Display query groups
    for query_group in query_groups.iter() {
        println!(
            "{}: {} sequences (edit cutoff: {})",
            query_group.barcode_type.as_str(),
            query_group.barcodes.len(),
            query_group.flank_k_cutoff.unwrap_or(0)
        );
    }
    println!();

    // Create demuxer
    let mut demuxer = Demuxer::new(
        0.05,  // alpha
        false, // verbose
        0.8,   // min_score
        0.1,   // min_score_diff
    );

    for query_group in query_groups.iter() {
        demuxer.add_query_group(query_group.clone());
    }

    // Load all reads into memory first
    println!("Loading reads into memory...");
    let mut reads: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
    let mut reader = parse_fastx_file(&fastq_path)?;

    while let Some(record) = reader.next() {
        let record = record?;
        reads.push((record.id().to_vec(), record.seq().to_vec()));
    }

    println!("Loaded {} reads", reads.len());
    println!();

    // Choose search strategy - enable all modes
    let mut search_strategy = SearchStrategy::new();
    search_strategy.enable(SearchMode::Adapter);
    search_strategy.enable(SearchMode::BarsAndFlanks);
    search_strategy.enable(SearchMode::JustBars);

    // Now run the profiling loop
    println!("Starting profiling loop...");
    let mut total_matches = 0;

    for iteration in 0..iterations {
        for (read_id, read_seq) in &reads {
            // Convert read_id from Vec<u8> to &str
            let read_id_str = std::str::from_utf8(read_id).unwrap_or("unknown");
            let matches = demuxer.demux(read_id_str, read_seq, &search_strategy);
            total_matches += matches.len();
        }

        if (iteration + 1) % (iterations / 10).max(1) == 0 {
            println!("Completed iteration {}/{}", iteration + 1, iterations);
        }
    }

    println!();
    println!("Profiling complete!");
    println!("Total demux calls: {}", reads.len() * iterations);
    println!("Total matches found: {}", total_matches);
    println!(
        "Average matches per read: {:.2}",
        total_matches as f64 / (reads.len() * iterations) as f64
    );

    Ok(())
}
