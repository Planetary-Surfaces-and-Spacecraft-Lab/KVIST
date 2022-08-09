#include "ist_spectra.hpp"
#include <tuple>

// Calculate zeros of Schrodinger scattering problem
ist_spectra_t ist_spectra(double Emin, double Emax, int numE, std::vector<double> xz, std::vector<double> u) {
    std::vector<double> Ez = allocate_linspace(Emin, Emax, numE);
    auto kz = vectorized_sqrt(Ez);

    // Intialize the following vectors with dummy values (1 & 0)
    std::vector M11vec(numE, 1.0);
    std::vector M12vec(numE, 1.0);
    std::vector M22vec(numE, 1.0);
    std::vector traceMvec(numE, 1.0);
    std::vector S(numE, 0.0);
    std::vector traceM(numE, 1.0);

    double M11prev = 1.0;

    // TODO: Can these be preallocated to a certain size? For now just let C++ standard library deal with that
    std::vector<double> M11signflipEz; 
    std::vector<double> M11signflipEzidx;
    std::vector<double> mu;
    std::vector<double> fmu;
    std::vector<double> Ej;
    std::vector<double> fEj;
    std::vector<double> Estar;
    std::vector<double> fEstar;
    std::vector<double> m;
    std::vector<double> a;    

    // This is created with a copy constructor - definitely doesn't need to be that
    // Also does direct pointer manipulation
    // TODO: Be Better
    class interval_search::M12 g(xz.data(), u.data(), xz.size());

    auto [M11, M12, M22, trMdiv2] = mmat_cpp_safe(Emin, xz, u);

    int N = mu.size();
    ist_spectra_t ret = {mu, fmu, Ej, fEj, Estar, fEstar, N, Ez, M12vec, M22vec, traceMvec, m, a};

    return ret;


}
