#include "math_util.h"
#include <math.h>
#include <chrono>
#include <iostream>
#include <numeric>


inlvector<double> allocate_linspace(double min, double max, const int n)
{
    auto t1 = std::chrono::high_resolution_clock::now();
    inlvector<double> ret(n);
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << ret[5] << std::endl;
    std::chrono::duration<double, std::micro> ms_double = t2 - t1;
    std::cout << "allocating time = " << ms_double.count() << "micros\n";

    t1 = std::chrono::high_resolution_clock::now();
    const double dx = (max - min) / (double(n) - 1.0);
    t2 = std::chrono::high_resolution_clock::now();
    ms_double = t2 - t1;
    std::cout << "scalar compute time = " << ms_double.count() << "micros\n";

    t1 = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < n; ++i)
    {
        ret[i] = min + (i * dx);
    }
    t2 = std::chrono::high_resolution_clock::now();
    ms_double = t2 - t1;
    std::cout << "computing time = " << ms_double.count() << "micros\n";
    return ret;
}

std::vector<double> vectorized_sqrt(std::vector<double> input)
{
    std::vector<double> output(input.size());
    int i = 0;

    for (int i = 0; i < input.size(); ++i)
    {
        output[i] = sqrt(input[i]);
    }

    return output;
}
