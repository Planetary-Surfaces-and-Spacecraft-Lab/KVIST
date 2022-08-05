#include "math_util.h"
#include <math.h>
#include <numeric>


std::vector<double> allocate_linspace(double min, double max, const int n)
{
    const double dx = (max - min) / (double(n) - 1.0);
    std::vector<double> ret(n);
    // ret.resize(n);
    
    for(int i=0; i<n; i++) {
        ret[i] = (min + (i*dx));
    }

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
