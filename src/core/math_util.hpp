#ifndef MATH_UTIL_HPP
#define MATH_UTIL_HPP
#include <vector>
#include <exception>
#include "interval_search.hpp"

double sign(double a);

double find_guess(interval_search::ZeroFunction& f, double lb, double ub);

// Implemented in math_util.cpp
std::vector<double> allocate_linspace(double min, double max, const int n);

// Implemented in math_util.cpp
std::vector<double> vectorized_sqrt(std::vector<double> input);

#endif /* !MATH_UTIL_HPP */