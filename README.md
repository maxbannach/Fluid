# Fluid

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3871709.svg)](https://doi.org/10.5281/zenodo.3871709)

This program implements a heuristic to find a shallow tree *T* for a given connected graph *G* such that *G* is contained in the closure of *T,* that is, a [treedepth decomposition](https://en.wikipedia.org/wiki/Tree-depth) of *G.*

The tool was developed by *Max Bannach, Sebastian Berndt, Martin Schuster* and *Marcel Wienöbst* as submission for the  heuristic track of [PACE 2020](https://pacechallenge.org/2020/). As such, the input format and I/O behavior is as specified by the PACE.

# Algorithm
The main idea is to use the recursive definition of treedepth in terms of minimal separators, that is, *td(G)* is the minimum *td(G\S)+|S|* over all minimal separators *S* of *G.* The heuristic tries to find a "good" separator by clustering the graph into communities using the [asynchronous fluid communities algorithm](https://arxiv.org/pdf/1703.09307.pdf). We then pick a separator that separates these communities as good as possible.

# Dependencies
The following open source [crates](https://crates.io) are used. They are automatically downloaded and compiled when the tool is build using *Cargo.*
- [rand](https://crates.io/crates/rand)
- [union-find](https://crates.io/crates/union-find)
- [signal-hook](https://crates.io/crates/signal-hook)

# Build
The heuristic is implemented in [Rust](https://www.rust-lang.org) and can simply be build using [Cargo](https://doc.rust-lang.org/cargo/getting-started/installation.html):

```
cargo build --release
```

# Run
After the build is completed, the tool can either be executed directly via

```
./target/release/fluid < <mygraph.gr>
```

or by using Cargo

```
cargo run --release < <mygraph.gr>
```

# Optil.io
For the PACE we submitted the static binary generated by `cargo build --release` to [optil.io](https://www.optil.io), that is, the file `./target/release/fluid`.
