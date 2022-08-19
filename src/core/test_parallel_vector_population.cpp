#include <vector>
#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>
#include <omp.h>

// compile with clang++ -O3 -Xclang -fopenmp -lomp test_parallel_vector_population.cpp -o vectest

// Some complicated x specific functions
double each_entry_compute(double x) {
    for(int i=0; i<100; ++i) {
        x = x + x/2.0;
    }
    return (x*x/6.0)+x;
}

using namespace std::chrono;

int main() {
    int veclength = 100000;

    // Preallocate some vectors
    std::vector<double> v1(veclength);
    std::vector<double> v2(veclength);
    std::vector<double> v3(veclength);
    std::vector<double> v4(veclength);

    // Method 1
    high_resolution_clock::time_point t1 = high_resolution_clock::now();
    for(int i=0; i<veclength; i++) {
        v1.push_back(each_entry_compute((double) i));
    }
    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double> >(t2 - t1);
    std::cout << "Method 1: " << time_span.count()*1e3 << " ms" << std::endl;

    // Method 2
    t1 = high_resolution_clock::now();
    for(int i=0; i<veclength; ++i) {
        v2[i] = each_entry_compute((double) i);
    }
    t2 = high_resolution_clock::now();
    time_span = duration_cast<duration<double> >(t2 - t1);
    std::cout << "Method 2: " << time_span.count()*1e3 << " ms" << std::endl;

    // Method 3
    // https://mac.r-project.org/openmp/
    t1 = high_resolution_clock::now();
    #pragma omp parallel for schedule(static)
    for(int i=0; i<veclength; ++i) {
        v3[i] = each_entry_compute((double) i);
    }
    t2 = high_resolution_clock::now();
    time_span = duration_cast<duration<double> >(t2 - t1);
    std::cout << "Method 3: " << time_span.count()*1e3 << " ms" << std::endl;

    return 0;
}