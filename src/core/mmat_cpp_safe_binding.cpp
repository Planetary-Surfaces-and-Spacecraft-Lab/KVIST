#include "mmat_cpp_safe_binding.hpp"
#include <exception>

/*
 * int mmat_automatic(const double E, const double* xz, const double* u, 
 *                    int numel, double* M, double* delx);
 */

std::tuple<double, double, double, double> mmat_cpp_safe(const double E, const std::vector<double> xz, const std::vector<double> u) {
    double M[4];
    double delx;
    int ret = mmat_automatic(E, xz.data(), u.data(), xz.size(), M, &delx);
    if (ret < 0) {
        throw std::runtime_error("mmat C failed to run");
    }

    double traceMdiv2 = (M[0]+M[2]) / 2.0;
    return std::make_tuple(M[0], M[1], M[3], traceMdiv2);
}
