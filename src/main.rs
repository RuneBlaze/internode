mod internode;
mod tree;
mod upgma;
use clap::{Parser, ArgEnum};
use internode::*;
use std::fs::{self, File};
use std::path::PathBuf;
use tracing::{info, warn};
use tree::{Mode, TreeCollection, UstarConfig};
use upgma::upgma_star;
use crate::tree::parse_newick;
use ndarray_npy::WriteNpyExt;

#[derive(Debug, ArgEnum, Clone, Copy)]
enum Preset {
    Vanilla,
    Abayes,
    HundredBootstrap,
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
    /// Range of the support threshold
    #[clap(short, long, parse(try_from_str = parse_bounds), default_value="0.0-1.0")]
    bounds: (f64, f64),
    /// Number of threads. Currently only useful for very large (2000+ genes and 50+ species) datasets.
    #[clap(short, long, default_value_t = 1usize)]
    threads: usize,
    /// Iteration method for distance imputation
    #[clap(long, default_value_t = String::from("uns"))]
    impute: String,
    /// Preset for the analysis
    #[clap(long, arg_enum)]
    preset: Option<Preset>,
    /// Experimental option, do not use yet.
    #[clap(long)]
    length_impute: bool,
    /// Only output the average distance matrix
    #[clap(long)]
    only_distances: bool,
}

fn parse_bounds(s: &str) -> Result<(f64, f64), String> {
    let mut parts = s.split('-');
    let a = parts
        .next()
        .ok_or("missing first bound")?
        .parse::<f64>()
        .expect("invalid first bound");
    let b = parts
        .next()
        .ok_or("missing second bound")?
        .parse::<f64>()
        .expect("invalid second bound");
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
    tracing_subscriber::fmt()
        .with_writer(std::io::stderr)
        .init();
    let mut config = args_to_config(&args);
    match args.preset {
        Some(Preset::Vanilla) => {
            config.mode = Mode::Internode;
        },
        Some(Preset::Abayes) => {
            config.mode = Mode::Support;
            config.lower_bound = 0.333;
            config.upper_bound = 1.0;
        },
        Some(Preset::HundredBootstrap) => {
            config.mode = Mode::Support;
            config.lower_bound = 0.0;
            config.upper_bound = 100.0;
        },
        None => {},
    }
    if let Some(preset) = args.preset {
        info!(?preset, "applied preset");
    }
    
    info!(
        "analysis started with mode {:?} using {} thread(s)",
        config.mode, args.threads
    );
    if config.mode == Mode::Support {
        info!(
            "support normalization scheme: linearly from [{}, {}] to [0, 1]",
            config.lower_bound, config.upper_bound
        );
    }
    let mut trees = TreeCollection::from_newick(args.input, &config).unwrap();
    info!(
        "read {} gene trees with {} taxa",
        trees.ngenes(),
        trees.ntaxa()
    );
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
    if args.only_distances {
        if let Some(out) = args.output {
            ustar.dm.write_npy(File::create(out)?)?;
        } else {
            // write to stdout
            ustar.dm.write_npy(std::io::stdout())?;
        }
        return Ok(()); // early return
    }
    let tree: String = if ustar.has_missing {
        info!(
            "found missing data, imputing missing distances in mode \"{}\"",
            args.impute
        );
        assert_eq!(args.impute, "uns");
        let upgma_tree = upgma_star(&ustar.dm, &ustar.mask)?;
        let impute_mode = if args.length_impute {
            Mode::NLength
        } else {
            Mode::Internode
        };
        impute_matrix(&mut ustar, &upgma_tree, impute_mode);
        let fastme_nni_tree_str =
            run_fastme(&trees.taxon_set, &ustar.dm, &FastMEConfig::new(true, false));
        let fastme_nni_tree = parse_newick(&mut trees.taxon_set, &fastme_nni_tree_str, &config);
        impute_matrix(&mut ustar, &fastme_nni_tree, impute_mode);
        ustar.raw_tree(&trees.taxon_set)
    } else {
        ustar.raw_tree(&trees.taxon_set)
    };
    if let Some(out) = args.output {
        fs::write(out, tree)?;
    } else {
        println!("{}", tree);
    }
    Ok(())
}
