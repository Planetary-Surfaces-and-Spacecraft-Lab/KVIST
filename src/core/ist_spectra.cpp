#include "ist_spectra.hpp"
#include <tuple>
#include <stdlib.h>

// Calculate zeros of Schrodinger scattering problem
ist_spectra_t ist_spectra(double Emin, double Emax, int numE, std::vector<double> xz, std::vector<double> u) {
    if(Emin >= Emax) {
        throw std::invalid_argument( "Emin >= Emax" );
    }
    std::vector<double> Ez = allocate_linspace(Emin, Emax, numE);
    auto kz = vectorized_sqrt(Ez);

    // Intialize the following vectors with dummy values (1 & 0)
    std::vector M11vec(numE, 1.0);
    std::vector M12vec(numE, 1.0);
    std::vector M22vec(numE, 1.0);
    std::vector traceMvec(numE, 1.0);
    std::vector S(numE, 0.0);
    std::vector traceM(numE, 1.0);

    double M11prev = 1.0;

    // TODO: Can these be preallocated to a certain size? For now just let C++ standard library deal with that
    std::vector<double> M11signflipEz; 
    std::vector<int> M11signflipEzidx;    

    // This is created with a copy constructor - definitely doesn't need to be that
    // Also does direct pointer manipulation
    // TODO: Be Better
    class interval_search::M12 g(xz.data(), u.data(), xz.size());

    auto [M11, M12, M22, trMdiv2] = mmat_cpp_safe(Emin, xz, u);

    int i = 2;
    for(auto E: Ez) {
        auto [M11, M12, M22, trMdiv2] = mmat_cpp_safe(E, xz, u);
        M11vec.push_back(M11);
        M12vec.push_back(M12);
        traceMvec.push_back(trMdiv2);
        
        auto delS = abs(sign(M11)-sign(M11prev))/2;
        if (delS > 0.0) {
            M11signflipEz.push_back(E);
            M11signflipEzidx.push_back(i);
        }
        S.push_back(S.back() + delS); // excluding nan case
        M11prev = M11;
        ++i;
    }
    // Ignore nan finding which is part of MATLAB version
    // The MATLAB version just ignores Nans which was done to make it work
    // This will do better (he says to himself)

    // Step I: Search for Estar which seperate degrees of freedom (N of these)
    int N = M11signflipEz.size();
    std::vector<double> Estar(N,0.0);
    std::vector<double> fEstar(N,0.0);
    double tolX = 1.0e-8; // Tolerance to use for bisection search
    auto M11f = interval_search::M11(xz.data(), u.data(), xz.size());

    // Do the bisection search from [Emin, M11signflipEz(0)] i.e. from Emin to first sign flip
    auto isres = interval_search::interval_search(Emin, M11signflipEz[0], M11f, tolX); 
    Estar.push_back(isres.x);
    fEstar.push_back(isres.feval);

    // Ideally we could speed this up somehow in parallel
    // std::vector may or may not be a good way to do it
    for(int i=1; i<N; i++) {
        auto lb = (M11signflipEz[i-1]+M11signflipEz[i])/2;
        // Last entry used Emax as upper bound
        auto ub = (i!=N-1) ? (M11signflipEz[i+1]+M11signflipEz[i])/2 : Emax; 
        isres = interval_search::interval_search(lb, ub, M11f, tolX);
        Estar[i] = isres.x;
        fEstar[i] = isres.feval;
    }

    // Step II: The Auxilary Spectrum
    std::vector<double> mu(N-1);
    std::vector<double> fmu(N-1);
    auto M12f = interval_search::M12(xz.data(), u.data(), xz.size());

    int midx = 0;
    for(i=1; i<N-1; i++) {
        auto lb = (i==2) ? Emin:Estar[i-1];
        auto ub = Estar[i];

        int signdiff = sign(M12f(lb)) != sign(M12f(ub));
        if(signdiff) {
            isres = interval_search::interval_search(lb, ub, M12f, tolX); 
            mu[midx] = isres.x;
            fmu[midx] = isres.feval;
            midx++;
        } else { // Possibly more than one zero crossing
            double guess = find_guess(M12f, lb, ub);
            if(guess == -1.0) {
                std::cerr << "No good guess. Skipping this loop\n";
                continue; // Skip this loop
            }

            if(sign(g(lb)) == sign(g(guess))) {
                std::cerr << "Bad guess, lb = " <<lb<<", guess = "<< guess << 
                             ", ub = "<< ub << std::endl;
            }

            // Left of guess
            isres = interval_search::interval_search(lb, guess, M12f, tolX); 
            mu[midx] = isres.x;
            fmu[midx] = isres.feval;
            midx++;

            // Right of guess
            isres = interval_search::interval_search(guess, ub, M12f, tolX); 
            mu[midx] = isres.x;
            fmu[midx] = isres.feval;
            midx++;
        }
    }

    // Step III: The Main Spectrum
    std::vector<double> Ej(2*N-1,0.0);
    std::vector<double> fEj(2*N-1,0.0);
    std::vector< std::pair<double, double>> EjfEj(2*N-1);
    auto trMp1f = interval_search::traceMdiv2_plus1(xz.data(), u.data(), xz.size());
    auto trMm1f = interval_search::traceMdiv2_minus1(xz.data(), u.data(), xz.size());
    isres = interval_search::interval_search(Ez[0], mu[0], trMp1f, tolX); 
    EjfEj[0] = std::make_pair(isres.x, isres.feval);
    isres = interval_search::interval_search(Ez[0], mu[0], trMm1f, tolX); 
    EjfEj[1] = std::make_pair(isres.x, isres.feval);

    for(int i=0; i<(N-2); i++) {
        isres = interval_search::interval_search(mu[i], mu[i+1], trMp1f, tolX); 
        EjfEj[2*i+1] = std::make_pair(isres.x, isres.feval);
        isres = interval_search::interval_search(mu[i], mu[i+1], trMm1f, tolX); 
        EjfEj[2*i+2] = std::make_pair(isres.x, isres.feval);
    }

    std::sort(EjfEj.begin(), EjfEj.end());
    double* M = (double*) malloc(4*sizeof(double));
    double* delx = (double*) malloc(1*sizeof(double));
    mmat_automatic(EjfEj.back().first,xz.data(),u.data(),xz.size(),M, delx);
    double lasttrM = (M[0]+M[1]);
    free(M);
    free(delx);

    // Do the last entry
    if(lasttrM < 0.0) {
        isres = interval_search::interval_search(mu.back(), Emax, trMp1f, tolX);
    } else {
        isres = interval_search::interval_search(mu.back(), Emax, trMm1f, tolX); 
    }
    EjfEj[2*N-2] = std::make_pair(isres.x, isres.feval);

    // Convert pairs to individual vectors
    for(int i=0; i<2*N-1; i++) {
        Ej[i] = EjfEj[i].first;
        fEj[i] = EjfEj[i].second;
    }

    // Calculate modulus and amplitudes of modes
    std::vector<double> m(N-1);
    std::vector<double> a(N-1);
    for(int j=0; j<N-1; j++) {
        m[j] = (Ej[2*j+2] - Ej[2*j+1]) / (Ej[2*j+2] - Ej[2*j]);
        a[j] = Ej[2*j+2] - Ej[2*j+1];
    }

    ist_spectra_t ret = {mu, fmu, Ej, fEj, Estar, fEstar, N, Ez, M12vec, M22vec, traceMvec, m, a};

    return ret;


}
