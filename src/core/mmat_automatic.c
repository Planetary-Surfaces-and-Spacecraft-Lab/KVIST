#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include <string.h>
#include "mmat_automatic.h"
// #include "mex.h" 

// Flip this to 1 to print extra information
#define DEBUG 0

#define DET_TOL 0.001


int mmat_automatic(const double E, const double* xz, const double* u, 
                   int N, double* M, double* delx){
    if(DEBUG) {
        printf("Starting mmat_automatic subroutine\n");
    }
    M[0] = 1.0; M[1] = 0.0; 
    M[2] = 0.0; M[3] = 1.0;

    if(DEBUG) {
        printf("Was able to access and set M\n");
    }
    double* Mtmp = (double*) malloc(4*sizeof(double));
    int real_kap_flag;
    *delx = fabs(xz[1]-xz[0]);

    if(DEBUG) {
        printf("Was able to access and set delx\n");
    }
    if(DEBUG) {
        printf("delx = %f\n",*delx);
        printf("N    = %i\n",N);
    }
    
    for(int i=1; i<N; i++) {
        double kap=0.0;
        double sig_plus_E = u[i]+E;
        
        // If signal big is negative, handle imaginary number        
        real_kap_flag = (sig_plus_E < 0.0) ? 0:1;
        kap = sqrt(fabs(sig_plus_E));
        
        // Domain error when calling sqrt above
        if(isnan(kap)) {
            free(Mtmp);
            return -1;
        }
        
        double T[4];
        //T[0] = 0.0; T[1]=1.0; T[2]=2.0; T[3]=0.0;
        
        
        if (real_kap_flag) {
            T[0] =  cos(kap*(*delx));     T[1] = (kap==0.0) ? 0.0:sin(kap*(*delx))/kap;
            T[2] = -kap*sin(kap*(*delx)); T[3] = cos(kap*(*delx));
        }
        else {
            T[0] =  cosh(kap*(*delx));     T[1] = (kap==0.0) ? 0.0:sinh(kap**(delx))/kap;
            T[2] =  kap*sinh(kap*(*delx)); T[3] = cosh(kap*(*delx));
        }
        
        if (DEBUG) {
            printf("kap flag = %i\n", real_kap_flag);
        }

        double detT = T[0]*T[3]-T[1]*T[2];
        
        if (DEBUG) {
            printf("i=%i\nkap = %f\nreal_kap_flag = %i\n",i,kap,real_kap_flag);
            printf("E = %f\nsig_plus_E = %f\n",E,sig_plus_E);
            printf("kap*(*delx) = %f\n",kap*(*delx)); 
            printf("T = %f %f\n    %f %f\n",T[0],T[1],T[2],T[3]);
            printf("det T = %f\ndet T - 1 = %f\n",detT,detT-1.0);
        }
        
        // det T should equal 1. If it doesn't to good numerical precision
        // that is a problem
        if (fabs(detT-1.0) >= DET_TOL || isnan(T[0])) {
            free(Mtmp);
            return -2;
        }
        

        // Calculate M=T*M with BLAS and stability scaling alg from:
        // https://arxiv.org/pdf/1507.00687.pdf
        double* Da = (double*) malloc(4*sizeof(double)); // diagonal matrix with max absolute value from each row
        Da[0] = (fabs(T[0]) > fabs(T[1])) ? fabs(T[0]):fabs(T[1]); Da[1] = 0.0;
        Da[2] = 0.0; Da[3] = (fabs(T[2]) > fabs(T[3])) ? fabs(T[2]):fabs(T[3]);
        double* Dainv = (double*) malloc(4*sizeof(double));
        double detDa = Da[0]*Da[3] - Da[1]*Da[2];
        Dainv[0] = Da[3]/detDa; Dainv[1] = 0.0;
        Dainv[2] = 0.0; Dainv[3] = Da[0]/detDa;
        double* Tp = (double*) malloc(4*sizeof(double)); // T prime;
        // Calculate Tp which is conditioned T
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 2, 2, 
                    1.0, Dainv, 2, T, 2, 0.0, Tp, 2);
        if(DEBUG) {
            printf("Da = %f %f\n       %f %f\n", 
                    Da[0],Da[1],Da[2],Da[3]);
            printf("Dainv = %f %f\n       %f %f\n", 
                    Dainv[0],Dainv[1],Dainv[2],Dainv[3]);
            printf("Tp = %f %f\n       %f %f\n", 
                    Tp[0],Tp[1],Tp[2],Tp[3]);
        }
        
        double* Db = (double*) malloc(4*sizeof(double)); // diagonal matrix with max absolute value from each row
        Db[0] = (fabs(M[0]) > fabs(M[1])) ? fabs(M[0]):fabs(M[1]); Db[1] = 0.0;
        Db[2] = 0.0; Db[3] = (fabs(M[2]) > fabs(M[3])) ? fabs(M[2]):fabs(M[3]);
        double* Dbinv = (double*) malloc(4*sizeof(double));
        double detDb = Db[0]*Db[3] - Db[1]*Db[2];
        Dbinv[0] = Db[3]/detDb; Dbinv[1] = 0.0;
        Dbinv[2] = 0.0; Dbinv[3] = Db[0]/detDb;
        double* Mp = (double*) malloc(4*sizeof(double)); // M prime;
        // Calculate Mp which is conditioned M
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 2, 2, 
                    1.0, M, 2, Dbinv, 2, 0.0, Mp, 2);
        
        if(DEBUG) {
            printf("M = %f %f\n    %f %f\n",M[0],M[1],M[2],M[3]);
            printf("Db = %f %f\n       %f %f\n", 
                    Db[0],Db[1],Db[2],Db[3]);
            printf("Mp = %f %f\n       %f %f\n", 
                    Mp[0],Mp[1],Mp[2],Mp[3]);
        }
                
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 2, 2, 
                    1.0, Tp, 2, Mp, 2, 0.0, Mtmp, 2);
        if(DEBUG) {
            printf("Mtmp(a) = %f %f\n          %f %f\n", 
                    Mtmp[0],Mtmp[1],Mtmp[2],Mtmp[3]);
            printf("Da = %f %f\n       %f %f\n", 
                    Da[0],Da[1],Da[2],Da[3]);
        }
        double* Mtmp2 = (double*) malloc(4*sizeof(double));
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 2, 2, 
                    1.0, Da, 2, Mtmp, 2, 0.0, Mtmp2, 2);
        if(DEBUG) {
            printf("Mtmp(b) = %f %f\n          %f %f\n", 
                    Mtmp2[0],Mtmp2[1],Mtmp2[2],Mtmp2[3]);
        }
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 2, 2, 
                    1.0, Mtmp2, 2, Db, 2, 0.0, Mtmp, 2);
        if(DEBUG) {
            printf("Mtmp(c) = %f %f\n          %f %f\n", 
                    Mtmp[0],Mtmp[1],Mtmp[2],Mtmp[3]);
        }
        cblas_dcopy(4,Mtmp,1,M,1);
        
        double detM = M[0]*M[3]-M[1]*M[2];

        if(isnan(M[0])) {
            printf("M is a nan\n");
            return -3;
        }
        
        if(DEBUG) {
            printf("Mtmp = %f %f\n       %f %f\n", 
                    Mtmp[0],Mtmp[1],Mtmp[2],Mtmp[3]);
            printf("M = %f %f\n    %f %f\n",M[0],M[1],M[2],M[3]);
            printf("det M = %f\n",detM);
            printf("\n");
        }
        
        free(Db);
        free(Dbinv);
        free(Mp);
        free(Da);
        free(Tp);
        free(Dainv);
        free(Mtmp2);
        
        // det(M) should be 1. Not being so could be a sign of numerical problems
        /*if (fabs(detM-1.0) >= DET_TOL*fabs(M[1])) {
            free(Mtmp);
            return -3;
        }*/
        

    }
    
    free(Mtmp);
    return 0;
    
}