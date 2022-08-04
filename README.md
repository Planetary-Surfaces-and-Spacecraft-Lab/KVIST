# **K**d**V** **I**nverse **S**cattering **T**ransform

## Description
Computes the Periodic Inverse Scattering Transform of the Kortweg-De Vries equation using the method outlined in Chapter 17 of [Osborne, A. (2010). Nonlinear Ocean Waves and the Inverse Scattering Transform. Elsevier Science](https://www.elsevier.com/books/nonlinear-ocean-waves-and-the-inverse-scattering-transform/osborne/978-0-12-528629-9). The algorithm consists of two parts:

1. Computing the 2x2 real scattering matrix $M$ of the discrete scattering potential $U(x)$ for a given energy level $E$.
2. Using a bisection search to compute the energy levels which are the zeros of $M_{12}$, $M_{11}$, and $\frac{1}{2} tr M$.

### Code Organization
- `src` subroutines to calculate $M$ and the zeros $E$. 
- `publications` directly reproduces figures associated with this code
- `dev_scripts` sample MATLAB scripts to call subroutines in `src`


## Installation
Thus far, this package has only been used on an Intel Macbook. In the future, more bindings and platforms are expected to be supported.

### Requirements
- MATLAB installation (with coder package). Must be in your `PATH`. 
- `make`
- C/C++ compiler (tested with Apple clang version 13.1.6)
- `BLAS` (for example: [OpenBLAS](https://www.openblas.net/))

### Instructions
In root directory run `make`. This will compile the necesary code and copy the entire library to `KVIST/lib`. In MATLAB, just add `KVIST/lib/` to your path and you're good to go! `make clean` will clean up the `src` folder. You make need to adjust the BLAS link and include facts in MATLAB coder. This can be done by editing `src/matlab/compile_mex.m`. Future versions of this will improve upon this.

The scripts in `publications` will search for a `lib` directory in order to run. 
