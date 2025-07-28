use benchmarks::simulations::cli::*;
use benchmarks::simulations::sim_data::*;
use clap::Parser;

fn main() {
    let cli = Cli::parse();
    println!("Output directory: {}", cli.output_dir);

    match cli.command {
        Commands::Simulate(args) => {
            println!("Running simulation with {} iterations", args.n);
            create_testdata(args.n, &args.output_dir, &args.barcode_file, args.rc_frac);
        }
    }
}
