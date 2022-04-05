function ret = M11f(E,xz,signal,lambda)
    ret = zeros(length(E),1);
    for i = 1:length(E)
        if length(E) > 1
            [M,~] = mmat_mex(E(i), xz, lambda*signal);
            ret(i) = M(1,1);
        else
            [M,~] = mmat_mex(E, xz, lambda*signal);
            ret = M(1,1);
        end
    end
end