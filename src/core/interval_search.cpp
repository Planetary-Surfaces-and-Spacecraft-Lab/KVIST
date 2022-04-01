#include "interval_search.hpp"
using namespace std;

// f(x) = x^2
class Squared : public ZeroFunction
{
public:
    double f(double x) {
        return x*x;
    }
};

// Given an interval where sign(f(a)) /= sign(f(b)), search for where f(x)=0
// Assumes a<b
IntervalSearchResult interval_search(double a, double b, ZeroFunction& f, double tolX) {
    double leftbnd = a;
    double rightbnd = b;
    double middle;
    double fhalf;

    // Make left half negative, right half positive
    double flip = (f.g(leftbnd)<=0.0)?1.0:-1.0;

    int iter = 0;
    int at_lftbnd = 1;
    int at_rhtbnd = 1;
    double delx = 1000.0; // High value until it gets set
    while(delx >= tolX) {
        middle = (leftbnd+rightbnd)/2.0;
        fhalf = f(middle);
        if (flip*fhalf < 0.0) { // Same sign as f(a)
            if (DEBUG){cout << "Changing left bound f(" <<middle<<") = "<<fhalf<<endl;}
            delx = middle - leftbnd;
            leftbnd = middle;
            at_lftbnd = 0;
        } 
        if (flip*fhalf > 0.0) { // Same sign as f(b)
            if(DEBUG){cout << "Changing right bound f(" <<middle<<") = "<<fhalf<<endl;}
            delx = rightbnd - middle;
            rightbnd = middle;
            at_rhtbnd = 0;
        }

        if (fhalf == 0.0) { // This is the solution
            struct IntervalSearchResult ret = {middle, fhalf, delx, iter};
            return ret;
        }

        iter++;
    }

    if(at_lftbnd) {
        if(DEBUG){cout << "Solutions is at left boundary" << endl;}
        middle = a;
        fhalf = f(a);
    }

    if(at_rhtbnd) {
        if(DEBUG){cout << "Solutions is at right boundary" << endl;}
        middle = b;
        fhalf = f(b);
    }

    struct IntervalSearchResult ret = {middle, fhalf, delx, iter};
    return ret;
}

// int main() {
//     Squared f;
//     IntervalSearchResult res1 = interval_search(0.0, 1.0, f, 1.0e-100);
//     cout<<"Results:\nx = "<<res1.x<<endl<<"f(x) = "<<res1.feval<<endl;
//     cout<<"Iterations: "<<res1.numiter<<endl<<endl;

//     IntervalSearchResult res2 = interval_search(-1.0, 0.0, f, 1.0e-100);
//     cout<<"Results:\nx = "<<res2.x<<endl<<"f(x) = "<<res2.feval<<endl;
//     cout<<"Iterations: "<<res2.numiter<<endl<<endl;

//     cout<<"Difference between xleft and xright = "<< res1.x-res2.x<<endl;
// }