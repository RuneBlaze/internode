mod internode;
use clap::Parser;
use internode::*;
use std::fs;
// use internode::tests::avian_tree;
use ndarray::prelude::*;
use std::path::PathBuf;

pub fn avian_tree() -> PathBuf {
    let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    d.push("resources/test");
    d.push("avian.tre");
    d
}

pub fn trivial_tree() -> PathBuf {
    let mut d = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    d.push("resources/test");
    d.push("trivial.tre");
    d
}

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
    /// Upper bound of the support values for wASTRID-s, the divider of the support
    #[clap(short = 'x', long, default_value_t = 1.0)]
    max_support: f64,
    /// Support value that denotes no information (support lower bound), for wASTRID-s
    #[clap(short, long, default_value_t = 0.0)]
    normalizer: f64,
    /// Number of threads. Currently only useful for very large (2000+ genes and 50+ species) datasets.
    #[clap(short, long, default_value_t = 1usize)]
    threads : usize,
}

fn args_to_config(args: &Args) -> UstarConfig {
    UstarConfig {
        max_support: args.max_support,
        normalizer: args.normalizer,
        mode: args.mode,
    }
}

fn main() {
    let args = Args::parse();
    let config = args_to_config(&args);
    let trees = TreeCollection::from_newick(args.input, &config).unwrap();
    let mut ustar = if args.threads == 1 {
        UstarState::from_tree_collection(&trees, &config)
    } else {
        rayon::ThreadPoolBuilder::new().num_threads(args.threads).build_global().unwrap();
        UstarState::from_tree_collection_par(&trees, &config, args.threads)
    };
    let tree = ustar.to_tree(&trees.taxon_set);
    if let Some(out) = args.output {
        fs::write(out, tree).unwrap();
    } else {
        println!("{}", tree);
    }
}