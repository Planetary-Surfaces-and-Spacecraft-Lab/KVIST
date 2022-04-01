function M11 = M11f(E,xz,signal,lambda)
    [M,~] = mmat_mex(E, xz, lambda*signal);
    M11 = M(1,1);
end