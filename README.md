# 3D-PDR
Code for modelling the chemistry of photodissociation regions in 1D and 3D

Click [here](http://itamos.readthedocs.io/) to visit the manual of the ITAMOS project and the 3D-PDR code.

## Quick installation 

### SUNDIALS

You will need to have `cmake` already installed in your system.

Install SUNDIALS using the following commands in a directory where you will generally have the 3D-PDR code.

```console
$ git clone https://github.com/LLNL/sundials.git sundials
```

Create and enter a build directory:

```console
$ mkdir sundials/build
$ cd sundials/build
```

Configure SUNDIALS using `cmake`:

```console
cmake -DCMAKE_INSTALL_PREFIX=/YOUR-HOMEPATH/3D-PDR/sundials ../
```

Build and install

```console
$ make
$ make install
```

Next, edit your shell configuration file (e.g. `~/.bashrc`) and add the following lines:

```console
export LD_LIBRARY_PATH=/home/USERNAME/3D-PDR/sundials/lib
export SUNDIALS_DIR=/home/USERNAME/3D-PDR/sundials
```

Then, reload your shell configuration:

```console
$ source ~/.bashrc
```

## PDR-studio

Navigate to the `PDR-studio/` directory inside your cloned `3D-PDR/` repository and run the following:

```bash
$ chmod 755 run.sh
$ ./run.sh
```

This will install all required dependencies. This step only runs once. Once the installation completes, you will be prompted with IP addresses. If your browser does not open automatically, Ctrl+Click on one of the displayed addresses.

You are now ready to run and analyse your models with PDR-studio.

> **Note:** The Chemical Analysis section requires an Ollama LLM model to generate AI summaries of your results. Instructions for installing and configuring Ollama are provided within that section.

## RAM estimator

`ram_estimate.py` (in the repository root, next to `params.dat`) estimates the peak RAM a 3D model will need **before you run it**, so you can size a job for a cluster/queue without guessing. It reads your `params.dat` (coolant list, HEALPix level) and `src/config.mk` (RAYTHEIA mode, NETWORK, CHEMANALYSIS, and the other compile flags that affect per-cell memory layout), fetches `NLEV`/`NTEMP` from each coolant's LAMDA file, and reports a per-component RAM breakdown for a grid of the resolution you specify.

It's a plain Python 3 script with no dependencies beyond the standard library (`gfortran` is used opportunistically to measure exact struct sizes — see `--no-probe` below).

### Usage

```console
$ python3 ram_estimate.py --res 256 --raytheia 1 --lmax 13
```

This reads `params.dat` in the current directory by default; pass `--params /path/to/params.dat` to point at a different model. `--res N` (required) sets the grid resolution (N³ cells). `--raytheia`, `--lmax`, and the chemistry network are normally picked up automatically from `src/config.mk`, but can all be overridden on the command line to explore "what if" scenarios (e.g. checking RAM at a higher HEALPix level, or with `-lmax` capping the number of coolant energy levels, before committing to a run).

Useful options:

| Option | Purpose |
|---|---|
| `--res N` | Grid resolution, N³ cells (required) |
| `--params FILE` | Path to `params.dat` (default: `params.dat`) |
| `--healpix L` | HEALPix level override (default: read from `params.dat`) |
| `--lmax N` | Cap energy levels per coolant, matching the `-lmax=N` runtime flag (default: no cap) |
| `--raytheia {0,1,2}` | Ray-tracing mode matching `config.mk RAYTHEIA=` (default: read from `config.mk`) |
| `--threads T` | OpenMP thread count, for the peak temporary (`evalpop`) estimate |
| `--chemanalysis {0,1}` | Override the `CHEMANALYSIS` flag (default: read from `config.mk`) |
| `--srcdir DIR` | Path to `src/`, used to locate `config.mk` and to measure exact struct sizes (default: `<params dir>/src`) |
| `--no-probe` | Skip compiling the struct-size probe with `gfortran`; use a rougher fallback instead |
| `--coolheat-fixed {0,1}` | Override detection of whether the `cooling`/`heating` allocation fix (see below) is present in `3DPDR.F90` |

Run `python3 ram_estimate.py --help` for the full list.

### What it reports

The tool prints, for the requested resolution:

- A per-coolant table of energy levels used (after any `-lmax` cap).
- A persistent-RAM breakdown across all components that scale with grid size: coolant level-population arrays, ray-tracing path arrays, chemical abundances, the Fortran `pdr_node`/`pdr_excit` struct and pointer-descriptor overhead (measured directly via `gfortran` when available), static LAMDA data, and — when relevant — the `CHEMANALYSIS` `temp_rate(nreac,pdr_ptot)` array, which is disproportionately large and explicitly flagged as "not recommended for 3D" in `config.mk`.
- The peak *temporary* RAM used during cooling-function evaluation (`evalpop`, scaled by `--threads`).
- A HEALPix-level comparison table, since ray-tracing memory scales as 12×4^L per cell and is usually the first thing to reconsider if a model doesn't fit in available RAM.

### Notes on accuracy

The estimate is derived directly from the Fortran source (`initialization.F90`, `allocations.F90`, `modules.F90`, `3DPDR.F90`), not from a fixed formula, so it stays correct as the code evolves. In particular:

- Per-cell struct overhead (array descriptors inside `pdr_node`/`pdr_excit`) is measured by compiling a small probe against your actual `src/modules.F90` with the same preprocessor flags as your build, rather than guessed.
- It auto-detects (and can be overridden) whether the `cooling(:)`/`heating(:)` allocation in `3DPDR.F90` uses the fixed, single correctly-sized `allocate` per cell, or the older per-index loop, which affects real memory usage.
- Expect the reported total to be within a few percent of observed usage; the remainder is typically OS/allocator overhead (heap alignment, malloc bookkeeping across the many small per-cell allocations) that isn't practical to model exactly, plus MPI buffers, Fortran runtime overhead, and CVODE solver workspace, none of which are included.

