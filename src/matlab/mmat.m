function [newM,k,delx] = mmat(E, xz, u)%#codegen
%#codegen
% for code generation, preinitialize the output variable
% data type, size, and complexity 
% generate an include in the C code
coder.cinclude('mmat_automatic.h');
M = zeros(4,1);
k = sqrt(complex(E));
% evaluate the C function
flag = int32(0);
N = int32(length(xz));
delx = 0.0;
coder.updateBuildInfo('addLinkFlags','-L /usr/local/opt/openblas/lib');
coder.updateBuildInfo('addLinkFlags','-lcblas');
coder.updateBuildInfo('addIncludePaths','/usr/local/opt/openblas/include');
flag = coder.ceval('mmat_automatic', E, coder.rref(xz), coder.rref(u),...
                    int32(N), coder.wref(M), coder.wref(delx)); 

if(flag < 0)
    error('Something went wrong, flag = %g\n',flag);
end
newM = zeros(2,2);
newM(1,1) = M(1);
newM(1,2) = M(2);
newM(2,1) = M(3);
newM(2,2) = M(4);

end