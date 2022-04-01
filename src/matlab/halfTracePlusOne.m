function Mtrdiv2 = halfTracePlusOne(E,xz,signal,lambda)
    [M,~] = mmat_mex(E, xz, lambda*signal);
    if(isnan((trace(M)/2) + 1))
        Mtrdiv2 = realmax;
    else
        Mtrdiv2 = (trace(M)/2) + 1;
    end
end