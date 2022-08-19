#include "math_util.hpp"
#include <math.h>
#include <numeric>
#include <stdexcept>
#include <cstdlib>
#include <random>
#include <iostream>

double find_guess(interval_search::ZeroFunction& f, double lb, double ub) {
    if(sign(f(lb)) != sign(f(ub))) {
        throw std::invalid_argument( "lower bound and upper bound need to be different arguments" );
    }
    double guess;

    // Try midway point
    if (sign(f(lb+ub)/2) != sign(f(lb))) {
        guess = (lb+ub)/2;
        return guess;
    }

    // Try quarterway point
    if (sign(f(lb+ub)/4) != sign(f(lb))) {
        guess = (lb+ub)/4;
        return guess;
    } 

    // Start randomly guessing
    int maxtry = 10000; // Limit number of guesses
    int it = 0;
    double rng_len = ub-lb;
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(lb, ub);
    guess = dis(gen) + lb;
    while(sign(f(guess) == sign(f(lb)))) {
        guess = dis(gen) + lb;
        it++;
        if(it >= maxtry) {
            std::cerr << "Could not find a good guess\n";
            guess = -1;
            return guess;
        }
    }

    return guess;
}

double sign(double a) {
    if(a>0.0) {
        return 1.0;
    } else if (a < 0.0){
        return -1.0;
    } else {
        return 0.0;
    }
}

std::vector<double> allocate_linspace(double min, double max, const int n)
{
    const double dx = (max - min) / (double(n) - 1.0);
    std::vector<double> ret;
    ret.resize(n);
    // ret.resize(n);
    
    for(int i=0; i<n; i++) {
        ret[i]=(min + (i*dx));
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
