# **K**d**V** **I**nverse **S**cattering **T**ransform

## Description

## Installation
Thus far, this package has only been used on an Intel Macbook. In the future, more bindings and platforms are expected to be supported.

### Requirements
- MATLAB installation (with coder package). Must be in your `PATH`
- `make`
- C/C++ compiler (tested with Apple clang version 13.1.6)
- `BLAS` (for example: [OpenBLAS](https://www.openblas.net/))

### Instructions
In root directory run `make`. This will compile the necesary code and copy the entire library to `KVIST/lib`. In MATLAB, just add 'KVIST/lib/` to your path and you're good to go! `make clean` will clean up the `src` folder.

The scripts in `publications` will search for a `lib` directory in order to run. 
