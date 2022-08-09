#include <vector>
#include "math_util.hpp"
#include "interval_search.hpp"
#include "mmat_cpp_safe_binding.hpp"

typedef struct ist_spectra_t
{
    std::vector<double> mu;
    std::vector<double> fmu;
    std::vector<double> Ej;
    std::vector<double> fEj;
    std::vector<double> Estar;
    std::vector<double> fEstar;
    int N;
    std::vector<double> Ez;
    std::vector<double> M12;
    std::vector<double> M22;
    std::vector<double> traceM;
    std::vector<double> m; // moduli of hyperelliptic modes
    std::vector<double> a; // amplitude of hyperelliptic modes
} ist_spectra_t;

ist_spectra_t ist_spectra(double Emin, double Emax, int numE, std::vector<double> xz, std::vector<double> u);