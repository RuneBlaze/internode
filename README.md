Internode
==================

[![shields.io](https://img.shields.io/badge/recommended_version-0.0.7_snapshot-blue?style=for-the-badge)](https://github.com/RuneBlaze/internode/releases/tag/v0.0.7-snapshot) [![shields.io](https://img.shields.io/badge/research_paper-wabi_2022-blue?style=for-the-badge)](https://drops.dagstuhl.de/opus/volltexte/2022/17042/)

Very fast and accurate species tree inference despite [ILS](https://en.wikipedia.org/wiki/Incomplete_lineage_sorting). Somewhat a successor to [ASTRID-2](https://github.com/pranjalv123/ASTRID) and a competitor to [ASTER](https://github.com/chaoszhang/ASTER).

## Quick Start

### Installing the binary

A single binary called `wastrid` can be downloaded from [our recommended release](https://github.com/RuneBlaze/internode/releases/tag/v0.0.7-snapshot), currently only containing binaries for x86_64 Linux and arm64 macOS. You can put it in `PATH` or simply in the directory that you wish to run it from.

### Preparing data

If you are familiar with [ASTRAL](https://github.com/smirarab/ASTRAL)/[ASTER](https://github.com/chaoszhang/ASTER), you can skip this section because we almost use the same input format.

`wastrid` accepts the common newline delimited Newick format for the input gene trees, where each line contains a gene tree in Newick format. We packed some example data from [S100](https://raw.githubusercontent.com/RuneBlaze/internode/main/resources/test/s100_k200.tre)[^1] containing 200 gene trees and 101 taxa.

In addition, the gene trees can be annotated with branch support that you trust (the example data given above is supplied with bootstrap support from 0 to 100) for better accuracy using wASTRID-s.

### Inferring a Species Tree

Say, you have prepared this `genes.tre` file containing all your gene trees separated by newlines. Now we have two differing scenarios:

 1. The gene trees have been annotated by support values that you trust (e.g., bootstrap support, [aBayes](https://academic.oup.com/sysbio/article/60/5/685/1644562?login=false) support). Moreover, you know the upper-bound and lower-bound of these support. Continuing with our example data of S100 mentioned above, the lower bound is 0 and the upper bound is 100.
 2. In the other scenario, the gene trees don't have support you trust.

In the first scenario, the command goes something like this:

```shell
wastrid -i genes.tre -b 0-100 -o output_stree.tre
```

where `output_stree.tre` is the output species tree path, `-b` specified our bounds for the support (lower bound 0, upper bound 100). For example, for aBayes support that can have lower bound 0.333 and upper bound 1, the bounds can be specified as `-b 0.333-1`.

In the second scenario (you don't have support/you don't trust the accuracy of such support), you really just want to run ASTRID (not weighted ASTRID), in which case just do

```shell
wastrid -i genes.tre --preset vanilla -o output_stree.tre
```

where `preset` can preconfigure flags for you. Other presets include:

 - `--preset abayes`, equivalent to `-m support -b 0.333-1`
 - `--preset hundred-bootstrap`, equivalent to `-m support -b 0-100`

After running the appropriate command, the output species tree topology is at `output_stree.tre`. Note that the branch lengths of the species tree are not biologically meaningful.

## Examples

Emulation of ASTRID: output is printed to STDOUT

```bash
wastrid -i gtrees.tre --preset vanilla
```

Weighted ASTRID by support on IQTree aBayes support: output is written to a file
specified by the path `species.tre` (`[]` denotes optional arguments)

```bash
wastrid -i gtrees.tre -b 0.333-1 [-m support] -o species.tre
```

## Compilation

See also prebuilt binaries (located in [Releases](https://github.com/RuneBlaze/internode/releases)).

`internode` is developed with Rust, so compiling it from scratch requires a proper installation of the [Rust toolchain](https://www.rust-lang.org/learn/get-started). However, the FastME bindings used internally make the compiling process a bit harder than simply `cargo build`.

The FastME bindings used internally were generated via `bindgen`, which
requires libclang. After a proper installation of `libclang` for `bindgen`, running the usual pipeline works for building the `wastrid` binary:

```bash
cargo build --release
```

For maximum speed that reduces cross-machine compatibility, the usual tricks apply:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release
```

## Notes
 - This implementation of ASTRID is faster than the original implementation (of ASTRID-2). That is, `wastrid --preset vanilla` is speed-wise a better ASTRID.
 - Missing data imputation is implemented (and automatically turned on), but alpha quality, using the original procedure of ASTRID.
 - ASTRID-multi (see also [DISCO](https://github.com/JSdoubleL/DISCO)) is still not implemented

## Acknowledgments

The code contains translated parts from [ASTRID-2](https://github.com/pranjalv123/ASTRID) and [TreeSwift](https://github.com/niemasd/TreeSwift). Due to ASTRID-2's GPLv2 and TreeSwift's GPLv3 license, this project is licensed as GPLv3.

FastME (residing in `third_party`) is derived from its
[original source code](https://gite.lirmm.fr/atgc/FastME/), and
of course, does not fall into the GPL licensing part of this project. I do not know a proper software license of FastME, so if this project
results in some violation please let me know.

`internode` and USTAR methods in general descend from the work of Liu and
collaborators. See [these](https://doi.org/10.1093/sysbio/syr027) [papers](https://doi.org/10.1186/1471-2164-16-S10-S3) [for](https://doi.org/10.1109/TCBB.2016.2604812) a start.

[^1]: From [ASTRAL-III](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2129-y), generated by Zhang and Mirarab.