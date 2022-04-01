addpath('../core/');

% Codegens and compiles mmat.m
E = 0.1;
xz = zeros(1,299); u = zeros(1,299);
eg_E = coder.typeof(E,[1],1);
max_size = 10000;
eg_xz = coder.typeof(xz,[1,max_size],1);
eg_u = coder.typeof(u,[1,max_size],1);
codegen mmat -args {eg_E, eg_xz, eg_u} ../core/mmat_automatic.c

% Command to mex cpp bindings for zero functions
mex ../core/BetterIntervalSearch_PlusOne.cpp ../core/mmat_automatic.c ../core/interval_search.cpp -lcblas -L/usr/local/opt/openblas/lib -I/usr/local/opt/openblas/include
mex ../core/BetterIntervalSearch_MinusOne.cpp ../core/mmat_automatic.c ../core/interval_search.cpp -lcblas -L/usr/local/opt/openblas/lib -I/usr/local/opt/openblas/include
mex ../core/BetterIntervalSearch_M12.cpp ../core/mmat_automatic.c ../core/interval_search.cpp -lcblas -L/usr/local/opt/openblas/lib -I/usr/local/opt/openblas/include
