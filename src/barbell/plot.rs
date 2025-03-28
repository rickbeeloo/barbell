
use plotly::{Plot, Scatter};




pub fn plot_edit_distances(all_edits: Vec<(usize, Vec<i32>)>) -> Result<(), Box<dyn std::error::Error>> {
    let mut plot = Plot::new();
    
    // Create a scatter plot for each flank's edits
    for (flank_idx, edits) in all_edits.iter() {
        let x: Vec<f64> = (0..edits.len()).map(|i| i as f64).collect();
        let y: Vec<f64> = edits.iter().map(|&e| e as f64).collect();
        
        let trace = Scatter::new(x, y)
            .name(format!("Flank {}", flank_idx))
            .mode(plotly::common::Mode::Lines);
        
        plot.add_trace(trace);
    }
    // Customize the layout
    plot.set_layout(plotly::Layout::new()
        .title(plotly::common::Title::from("Edit Distances Across Read Positions"))
        .x_axis(plotly::layout::Axis::new().title(plotly::common::Title::from("Position")))
        .y_axis(plotly::layout::Axis::new().title(plotly::common::Title::from("Edit Distance"))));
    
    plot.show();
    
    Ok(())
}