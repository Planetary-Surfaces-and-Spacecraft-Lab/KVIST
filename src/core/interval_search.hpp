#ifndef INTERVAL_SEARCH_HPP
#define INTERVAL_SEARCH_HPP
#include <stdlib.h>
#include <iostream>
extern "C" {
    #include "mmat_automatic.h"
}

#define DEBUG 0

namespace interval_search {
// This would be best as a template to allow f(x) other than double |-> double
class ZeroFunction {
public:
    virtual double f(double a)=0;
    virtual double g(double a)=0;

    // Allows () operator to return f(a)
    // This pretty similar to how functions are implemented in Python I think
    double operator()(double a) {
        return f(a);
    }
};

typedef struct IntervalSearchResult
{
    double x;
    double feval;
    double delx;
    int numiter;
}IntervalSearchResult;


// (tr M / 2) + 1
class traceMdiv2_plus1 : public ZeroFunction
{
private:
    double* xz;
    double* u;
    int n;
public:

    // Copy data into instance of ZeroFunction
    traceMdiv2_plus1(const double* xz_arg, const double* u_arg, int numel) {
        xz = (double *) malloc(numel*sizeof(double));
        u = (double *) malloc(numel*sizeof(double));
        for(int i=0; i<numel;i++) {
            xz[i] = xz_arg[i];
            u[i] = u_arg[i];
        }

        n = numel;
    }
    ~traceMdiv2_plus1() {
        free((void *)xz);
        free((void *)u);
    }
    double f(double E) {
        double* M = (double*) malloc(4*sizeof(double));
        double* delx = (double*) malloc(1*sizeof(double));
        mmat_automatic(E,xz,u,n,M, delx);
        double trMdiv2 = (M[0]+M[3])/2.0;
        free(M);
        free(delx);
        return trMdiv2+1.0;
    }

    double g(double E) {
        double* M = (double*) malloc(4*sizeof(double));
        double* delx = (double*) malloc(1*sizeof(double));
        mmat_automatic(E,xz,u,n,M, delx);
        double trMdiv2 = (M[0]+M[3])/2.0;
        free(M);
        free(delx);
        return trMdiv2;
    }
    
};

// (tr M / 2) - 1
class traceMdiv2_minus1 : public ZeroFunction
{
private:
    double* xz;
    double* u;
    int n;
public:

    // Copy data into instance of ZeroFunction
    traceMdiv2_minus1(const double* xz_arg, const double* u_arg, int numel) {
        xz = (double*) malloc(numel*sizeof(double));
        u = (double*) malloc(numel*sizeof(double));
        for(int i=0; i<numel;i++) {
            xz[i] = xz_arg[i];
            u[i] = u_arg[i];
        }
        n = numel;
    }
    ~traceMdiv2_minus1() {
        free(xz);
        free(u);
    }
    double f(double E) {
        double* M = (double*) malloc(4*sizeof(double));
        double* delx = (double*) malloc(1*sizeof(double));
        mmat_automatic(E,xz,u,n,M, delx);
        double trMdiv2 = (M[0]+M[3])/2.0;
        free(M);
        free(delx);
        return trMdiv2-1.0;
    }

    double g(double E) {
        double* M = (double*) malloc(4*sizeof(double));
        double* delx = (double*) malloc(1*sizeof(double));
        mmat_automatic(E,xz,u,n,M, delx);
        double trMdiv2 = (M[0]+M[3])/2.0;
        free(M);
        free(delx);
        return trMdiv2;
    }
};


// M_{12}
class M12 : public ZeroFunction
{
private:
    double* xz;
    double* u;
    int n;
public:

    // Copy data into instance of ZeroFunction
    M12(const double* xz_arg, const double* u_arg, int numel) {
        xz = (double *) malloc(numel*sizeof(double));
        u = (double *) malloc(numel*sizeof(double));
        for(int i=0; i<numel;i++) {
            xz[i] = xz_arg[i];
            u[i] = u_arg[i];
        }

        n = numel;
    }
    ~M12() {
        free((void *)xz);
        free((void *)u);
    }
    double f(double E) {
        double* M = (double*) malloc(4*sizeof(double));
        double* delx = (double*) malloc(1*sizeof(double));
        mmat_automatic(E,xz,u,n,M, delx);
        double M12 = M[1];
        free(M);
        free(delx);
        return M12;
    }

    double g(double E) {
        double* M = (double*) malloc(4*sizeof(double));
        double* delx = (double*) malloc(1*sizeof(double));
        mmat_automatic(E,xz,u,n,M, delx);
        double M12 = M[1];
        free(M);
        free(delx);
        return M12;;
    }
    
};

// M_{11}
class M11 : public ZeroFunction
{
private:
    double* xz;
    double* u;
    int n;
public:

    // Copy data into instance of ZeroFunction
    M11(const double* xz_arg, const double* u_arg, int numel) {
        xz = (double *) malloc(numel*sizeof(double));
        u = (double *) malloc(numel*sizeof(double));
        for(int i=0; i<numel;i++) {
            xz[i] = xz_arg[i];
            u[i] = u_arg[i];
        }

        n = numel;
    }
    ~M11() {
        free((void *)xz);
        free((void *)u);
    }
    double f(double E) {
        double* M = (double*) malloc(4*sizeof(double));
        double* delx = (double*) malloc(1*sizeof(double));
        mmat_automatic(E,xz,u,n,M, delx);
        double M11 = M[0];
        free(M);
        free(delx);
        return M11;
    }

    double g(double E) {
        double* M = (double*) malloc(4*sizeof(double));
        double* delx = (double*) malloc(1*sizeof(double));
        mmat_automatic(E,xz,u,n,M, delx);
        double M11 = M[0];
        free(M);
        free(delx);
        return M11;;
    }
    
};


IntervalSearchResult interval_search(double a, double b, ZeroFunction& f, double tolX);

};

#endif /* !INTERVAL_SEARCH_HPP */
