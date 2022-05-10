internode
==================

Fast implementation of ASTRID-like methods, alpha quality.

## Usage

The compiled binary name is `wastrid` (weighted ASTRID).

```bash
wastrid -i <input> -m mode [-o <output>]
```

where `<input>` is newline delimited single-copy gene trees in
newick format. `<output>` is the full path to the output species
tree, when unspecified, goes to `STDOUT`.

`mode` defaults to `support`, which sums up the internode support. Set it to `internode` for an emulation of the original ASTRID's behavior -- simply summing up `1` for each branch. If the support values are upper-bounded by somthing greater than 1 (say, if using bootstrap support in the values between 0 and 100), then the `-x` option needs to be specified to be the upper-bound (`-x 100`).

Likewise, if some support value should mean "no information" (AFAIK, aBayes from IQTree),
it should be specified as such via `-n`, for example `-n 0.333`. This way the support values are normalized from `[0.333, 1]` to `[0, 1]`. These options names are borrowed from weighted ASTRAL. To summarize, for support-weighted use-cases,

 - `wastrid -i $genes -n 0.333 -o $output` (aBayes)
 - `wastrid -i $genes -x 100 -o $output` (bootstrap support with support value upper-bound `100`)
 - `wastrid -i $genes -o $output` (any support that  has range `[0, 1]`, for example, the default support of FastTree)

Currently the support value normalization is different from what is used in weighted ASTRAL (as of `astral-hybrid` v1.4.2.3). The upper-bound (`-x`) is first divided then the lower bound (`-n`) subtracted, but for the above three types of support the flags are consistent across wASTRAL and this project.

For determining which support to use, aBayes support *seems* like the most accurate support to use (for weighted ASTRAL too AFAIK).

Note that missing data in the internode distance is not handled properly: each pair of taxa must appear in some gene tree.

## Examples

Emulation of ASTRID: output is printed to STDOUT

```bash
wastrid -i gtrees.tre -m internode
```

Weighted ASTRID by support on IQTree aBayes support: output is written to a file
specified by the full-path `species.tre`

```bash
wastrid -i gtrees.tre -n 0.333 [-m support] -o species.tre
```

## Building

Prebuilt binaries for x86_64 Linux (using the `x86_64-unknown-linux-musl` target) and M1 Mac (M1 Mac binaries not tested) are available in the Github releases.

`internode` is developed with Rust, so compiling it from scratch requires a proper installation of the [Rust toolchain](https://www.rust-lang.org/learn/get-started), but the FastME bindings used internally make the compiling process a bit harder than simply `cargo build`.

The FastME bindings used internally were generated via `bindgen`, which
requires libclang. After a proper installation of `libclang` for `bindgen`, running the usual pipeline works for building the `wastrid` binary:

```bash
cargo build --release
```

For maximum speed that reduces cross-machine compatibility, the usual tricks apply:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

## Missing features

 - Iteration ("missing data imputation") via UPGMA*. For each pair of taxa, it must appear simutaneously in some gene tree.
 - ASTRID-multi (see also [DISCO](https://github.com/JSdoubleL/DISCO))

## How fast?

```
░▒▓ ~/scratch/rigor/model.200.2000000.0.0000001/01 hyperfine --warmup 1 "ASTRID -s -i gtrees.tre.k1000" "wastrid -i gtrees.tre.k1000 -m internode -o gtrees.tre.k1000.internode"
Benchmark 1: ASTRID -s -i gtrees.tre.k1000
  Time (mean ± σ):      2.099 s ±  0.034 s    [User: 1.924 s, System: 0.174 s]
  Range (min … max):    2.061 s …  2.155 s    10 runs

Benchmark 2: wastrid -i gtrees.tre.k1000 -m internode -o gtrees.tre.k1000.internode
  Time (mean ± σ):     431.4 ms ±   3.1 ms    [User: 410.3 ms, System: 20.6 ms]
  Range (min … max):   428.3 ms … 437.7 ms    10 runs

Summary
  'wastrid -i gtrees.tre.k1000 -m internode -o gtrees.tre.k1000.internode' ran
    4.87 ± 0.09 times faster than 'ASTRID -s -i gtrees.tre.k1000'
```

## Acknowledgements

The code contains translated parts from [ASTRID-2](https://github.com/pranjalv123/ASTRID) and [TreeSwift](https://github.com/niemasd/TreeSwift). Due to ASTRID-2's GPLv2 and TreeSwift's GPLv3 license, this project is licensed as GPLv3.

FastME (residing in `third_party`) is derived from its
[original source code](https://gite.lirmm.fr/atgc/FastME/), and
of course, does not fall into the GPL licensing part of this project. I do not know a proper software license of FastME, so if this project
results in some violation please let me know.