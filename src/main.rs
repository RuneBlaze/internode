mod internode;
mod tree;
mod upgma;
use clap::Parser;
use internode::*;
use upgma::upgma_star;
use std::fs;
use std::path::PathBuf;
use tracing::{info, warn};
use tree::{Mode, TreeCollection, UstarConfig};

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Path to the input newline delimited gene trees
    #[clap(short, long)]
    input: PathBuf,
    /// Path to the output species tree topology
    #[clap(short, long)]
    output: Option<PathBuf>,
    /// Analysis mode
    #[clap(short, long, arg_enum, default_value_t = Mode::Support)]
    mode: Mode,
    /// Range of the support threshold
    #[clap(short, long, parse(try_from_str = parse_bounds), default_value="0.0 1.0")]
    bounds : (f64, f64),
    /// Number of threads. Currently only useful for very large (2000+ genes and 50+ species) datasets.
    #[clap(short, long, default_value_t = 1usize)]
    threads: usize,
    /// Iteration method for distance imputation
    #[clap(long, default_value_t = String::from("uns"))]
    impute: String,
}

fn parse_bounds(s: &str) -> Result<(f64, f64), String> {
    let mut parts = s.split_ascii_whitespace();
    let a = parts.next().ok_or("missing first bound")?.parse::<f64>().expect("invalid first bound");
    let b = parts.next().ok_or("missing second bound")?.parse::<f64>().expect("invalid second bound");
    if a > b {
        Err("first bound must be less than second bound".to_string())
    } else {
        Ok((a, b))
    }
}

fn args_to_config(args: &Args) -> UstarConfig {
    UstarConfig {
        upper_bound: args.bounds.1,
        lower_bound: args.bounds.0,
        mode: args.mode,
    }
}

fn main() -> anyhow::Result<()> {
    let args = Args::parse();
    let config = args_to_config(&args);
    tracing_subscriber::fmt::init();
    info!("analysis started with mode {:?} using {} threads", config.mode, args.threads);
    if config.mode == Mode::Support {
        info!("support normalization scheme: linearly from [{}, {}] to [0, 1]", config.lower_bound, config.upper_bound);
    }
    let trees = TreeCollection::from_newick(args.input, &config).unwrap();
    info!("read {} gene trees with {} taxa", trees.ngenes(), trees.ntaxa());
    let mut ustar = if args.threads == 1 {
        UstarState::from_tree_collection(&trees, &config)
    } else {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .unwrap();
        UstarState::from_tree_collection_par(&trees, &config, args.threads)
    };
    ustar.flatten();
    let tree : String = if ustar.has_missing {
        info!("Found missing data, imputing missing distances in mode \"{}\"", args.impute);
        let upgma_tree = upgma_star(&ustar.dm, &ustar.mask)?;
        impute_matrix(&mut ustar, &upgma_tree, Mode::Internode);
        ustar.raw_tree(&trees.taxon_set)
    } else {
        ustar.raw_tree(&trees.taxon_set)
    };
    if let Some(out) = args.output {
        fs::write(out, tree).unwrap();
    } else {
        println!("{}", tree);
    }
    Ok(())
}
