#include "mex.hpp"
#include "mexAdapter.hpp"
#include "interval_search.hpp"

using namespace matlab::data;
using matlab::mex::ArgumentList;

// Inputs: xz, u, Emin, Emax, tolX
class MexFunction : public matlab::mex::Function {
public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        checkArguments(outputs, inputs);
        matlab::data::ArrayFactory factory;
        size_t n = inputs[0].getNumberOfElements();
        TypedArray<double> xz_mtlb = std::move(inputs[0]);
        TypedArray<double> u_mtlb = std::move(inputs[1]);
        buffer_ptr_t<double> xz_buff_ptr = xz_mtlb.release();
        buffer_ptr_t<double> u_buff_ptr = u_mtlb.release();
        double* xz = xz_buff_ptr.get();
        double* u = u_buff_ptr.get();
        M12 f(xz, u, n);
        double Emin = inputs[2][0];
        double Emax = inputs[3][0];
        double tolX = inputs[4][0];
        
        IntervalSearchResult res = interval_search(Emin, Emax, f, tolX);
        
        outputs[0] = factory.createScalar<double>(res.x);
        outputs[1] = factory.createScalar<double>(res.feval);
        outputs[2] = factory.createScalar<double>(res.delx);
        outputs[3] = factory.createScalar<int>(res.numiter);
        return;
        
    }

    void checkArguments(ArgumentList outputs, ArgumentList inputs) {
        // TODO: Add sanity checks here
        return; 
    }
};