#include <iostream>
#include "ist_spectra.hpp"
#include "math_util.hpp"

int main() {
    auto xz = allocate_linspace(0.0, 10.0, 1001);
    auto u  = allocate_linspace(1.0, 1.0, xz.size());

    double Emin = 0.0;
    double Emax = 0.0;
    int numE = 1001;
    auto istspc = ist_spectra(Emin, Emax, numE, xz, u);


    return 0;
}
