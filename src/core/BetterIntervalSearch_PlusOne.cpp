#include "mex.hpp"
#include "mexAdapter.hpp"
#include "interval_search.hpp"
#include <iostream>

using namespace interval_search;
using namespace matlab::data;
using matlab::mex::ArgumentList;


// Inputs: xz, u, Emin, Emax, tolX
class MexFunction : public matlab::mex::Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        checkArguments(outputs, inputs);
        matlab::data::ArrayFactory factory;
        size_t n = inputs[0].getNumberOfElements();
        //std::cout << "Made a MATLAB factory\n";
        TypedArray<double> xz_mtlb = std::move(inputs[0]);
        //std::cout << "Made xz typedarray factory\n";
        TypedArray<double> u_mtlb = std::move(inputs[1]);
        //std::cout << "Made u typedarray factory\n";
        buffer_ptr_t<double> xz_buff_ptr = xz_mtlb.release();
        //std::cout << "Made xz buff_ptr factory\n";
        buffer_ptr_t<double> u_buff_ptr = u_mtlb.release();
        //std::cout << "Made u buff_ptr factory\n";
        double* xz = xz_buff_ptr.get();
        double* u = u_buff_ptr.get();
//         std::cout << "xz[0] = "<<xz[0]<<"\n";
//         std::cout << "u[0] = " << u[0]<<"\n";
//         std::cout << "n = " <<n<<"\n";
//         std::cout << "Got xz and u pointer\nCreating f object\n";
        traceMdiv2_plus1 f(xz, u, n);
        //std::cout << "Created f object\n";
        double Emin = inputs[2][0];
        //std::cout << "Emin = "<<Emin<<std::endl;
        double Emax = inputs[3][0];
        //std::cout << "Emax = "<<Emax<<std::endl;
        double tolX = inputs[4][0];
        //std::cout << "tolX = "<<tolX<<std::endl;
        
        //std::cout << "Calling interval search\n";
        IntervalSearchResult res = interval_search::interval_search(Emin, Emax, f, tolX);
        //std::cout << "Finished interval search with " <<res.numiter<<" iterations\n";
        
        outputs[0] = factory.createScalar<double>(res.x);
        outputs[1] = factory.createScalar<double>(res.feval);
        outputs[2] = factory.createScalar<double>(res.delx);
        outputs[3] = factory.createScalar<int>(res.numiter);
        
        //std::cout << "Set outputs\n";
        return;
        
    }

    void checkArguments(ArgumentList outputs, ArgumentList inputs) {
        // TODO: Add sanity checks here
        return; 
    }
};