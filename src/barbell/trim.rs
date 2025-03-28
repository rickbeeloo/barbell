
use crate::barbell::filter::AnnotationLine;




pub fn trim_matches(filtered_match_file: &str, read_file: &str, output_file: &str) {

    
    
    // let mut reader = csv::ReaderBuilder::new()
    //     .delimiter(b'\t')
    //     .from_path(filtered_match_file).expect("Failed to open filtered match file");


    // for result in reader.deserialize() {
    //     let record: AnnotationLine = result?;
        
    //     if let Some(read_id) = &current_read_id {
    //         if *read_id != record.read {
    //             // Process previous group
    //             if check_filter_pass(&current_group, &filters) {
    //                 kept_reads += 1;
    //                 for annotation in &current_group {
    //                     writer.write_record(&annotation.to_record())?;
    //                 }
    //             }
    //             current_group.clear();
    //             current_read_id = Some(record.read.clone());
    //         }
    //     } else {
    //         current_read_id = Some(record.read.clone());
    //     }
        
    //     current_group.push(record);
    //     total_reads += 1;
    //     progress_bar.set_message(format!("{}", total_reads));
    // }
    
    // // Process the last group
    // if !current_group.is_empty() && check_filter_pass(&current_group, &filters) {
    //     kept_reads += 1;
    //     for annotation in &current_group {
    //         writer.write_record(&annotation.to_record())?;
    //     }
    // }



  
}