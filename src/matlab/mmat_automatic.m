function [M,k] = mmat_automatic(E, xz, u)
%mmat Calculates spectral matrix for a given k
%   From equations 3.2 - 3.14 in 
%    Nonlinear Fourier Analysis for the 
%    Infinite-Interval Korteweg-de Vries 
%    Section 17.6 of Nonlinear Ocean Waves & Inverse Scattering Transform
%    "Automatic Numerical IST Algorithm"

k = sqrt(E);
N = length(xz);
delx = xz(2)-xz(1);

M = eye(2);

for n = 1:N
    kap = sqrt(u(n)+E);
    kpp = sqrt(abs(kap^2));
    if kap^2 > 0
        T = [cos(kpp*delx)      sin(kpp*delx)/kpp; ...
            -kpp*sin(kpp*delx)  cos(kpp*delx)];
    else
        T = [cosh(kpp*delx)      sinh(kpp*delx)/kpp; ...
            -kpp*sinh(kpp*delx)  cosh(kpp*delx)];
    end
        
    M = T*M;
end

end
