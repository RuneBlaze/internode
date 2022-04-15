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

`mode` defaults to `support`, which sums up the internode support. Set it to `internode` for an emulation of the original ASTRID's behavior. If your support values are upper-bounded by somthing greater than 1 (say, if you have bootstrap support in the values between 0 and 100), then the `-x` option needs to be specified to be the upper-bound (`-x 100`).

## Examples

Emulation of ASTRID: output is printed to STDOUT

```bash
wastrid -i gtrees.tre -m internode
```

Weighted ASTRID by support: output is written to a file
specified by the full-path `species.tre`

```bash
wastrid -i gtrees.tre -m support -o species.tre
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

This code contains translated parts from [ASTRID-2](https://github.com/pranjalv123/ASTRID) and [TreeSwift](https://github.com/niemasd/TreeSwift).