#include <iostream>
#include "ist_spectra.hpp"
#include "math_util.hpp"
#include <chrono>

using namespace std::chrono;


int main() {
    auto xz = allocate_linspace(0.0, 10.0, 1001);
    auto u  = allocate_linspace(1.0, 1.0, xz.size());

    double Emin = 0.0;
    double Emax = 100.0;
    int numE = 11;

    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    auto istspc = ist_spectra(Emin, Emax, numE, xz, u);
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double> >(t2 - t1);
    std::cout << "Time taken: " << time_span.count()*1e3 << " ms" << std::endl;

    return 0;
}
