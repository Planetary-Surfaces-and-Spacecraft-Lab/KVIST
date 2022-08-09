#include <tuple>
#include <vector>
#include "mmat_automatic.h"

// Returns M11, M12, and trace M / 2
std::tuple<double, double, double, double> mmat_cpp_safe(const double E, std::vector<double> xz, std::vector<double> u);
