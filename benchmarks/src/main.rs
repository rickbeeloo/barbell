use benchmarks::compare::cli::*;
use benchmarks::compare::compare::run_all_tools;
use benchmarks::simulations::cli::*;
use benchmarks::simulations::sim_data::*;
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    /// Output directory for results
    #[arg(short, long, default_value = "output")]
    pub output_dir: String,

    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Run simulations
    #[command(name = "sim")]
    Simulate(SimulationArgs),

    /// Run comparison
    #[command(name = "compare")]
    Compare(CompareArgs),
}

fn main() {
    let cli = Cli::parse();
    println!("Output directory: {}", cli.output_dir);

    match cli.command {
        Commands::Simulate(args) => {
            println!("Running simulation with {} iterations", args.n);
            create_testdata(args.n, &args.output_dir, &args.barcode_file, args.rc_frac);
        }
        Commands::Compare(args) => {
            println!("Running comparison with {} threads", args.threads);
            run_all_tools(
                &args.fastq_file,
                &args.output_dir,
                &args.bar_file,
                args.threads,
                Some(args.extra_file),
                &args.dorado_exec_path,
                &args.barbell_exec_path,
                &args.flexiplex_exec_path,
            );
        }
    }
}
