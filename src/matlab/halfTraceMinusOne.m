function Mtrdiv2 = halfTraceMinusOne(E,xz,signal,lambda)
    N = length(E);
    Mtrdiv2 = zeros(N,1);
    for i =1:N
        if N > 1
            [M,~] = mmat_mex(E(i), xz, lambda*signal);
        else
            [M,~] = mmat_mex(E, xz, lambda*signal);
        end
        if(isnan((trace(M)/2) - 1))
            Mtrdiv2(i) = realmax;
        else
            Mtrdiv2(i) = (trace(M)/2) - 1;
        end
    end
    
end