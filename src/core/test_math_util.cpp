#include <iostream>
#include <chrono>
#include <thread>
#include "math_util.hpp"

int main() {
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    auto t1 = high_resolution_clock::now();
    std::vector<double> v =  allocate_linspace(0.0, 9.0, 1000);
    auto t2 = high_resolution_clock::now();
    std::cout << "v (allocate):  ";
    std::cout << std::endl << "v[4] = " << v[4] << std::endl;
    // std::cout << "length = " << v.size();
    // for(auto a: v)
    //     std::cout << a << " ";
    std::cout << std::endl;
    auto ms_int = duration_cast<std::chrono::microseconds>(t2 - t1);

    /* Getting number of milliseconds as a double. */
    duration<double, std::micro> ms_double = t2 - t1;

    std::cout << ms_double.count() << " micro seconds\n";
    std::vector<double> vv(v.size());//) = v.stdvec();

    t1 = high_resolution_clock::now();
    std::vector<double> w = vectorized_sqrt(vv);
    t2 = high_resolution_clock::now();
    std::cout << "w (vectorized sqrt):  ";
    // for(auto a: w)
    //     std::cout << a << " ";
    std::cout << std::endl;
       /* Getting number of milliseconds as an integer. */
    ms_int = duration_cast<std::chrono::microseconds>(t2 - t1);

    /* Getting number of milliseconds as a double. */
    ms_double = t2 - t1;

    std::cout << ms_double.count() << " micro seconds\n";

    
    return 0;
}