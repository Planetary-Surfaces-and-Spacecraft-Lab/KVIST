clear all; close all; clc;

L = 250;
num_points = 5001;
x = linspace(-L, L, num_points);
delx = x(2)-x(1);
a = 100;
eta0 = 0.5;
%signal = zeros(size(x));
%signal(25:75) = 0.1;
%signal = sin(2*2*pi*x)+sin(3*2*pi*x)+cos(5*2*pi*x);
signal = zeros(size(x));
signal(abs(x)<a) = eta0;

fh = four(signal);
%kz = 0:2*pi/(delx*num_points):2*pi/delx-2*pi/(delx*num_points);
kz = linspace(0, 0.2, 5000);
h = 10;
lambda = 1;%abs(eta0/h); % Does this match Al example?
[S,b,Fc,M11,M12,traceM,Ez] = condst(kz, x, signal, lambda);
kinf = 1e18;
minfpos11 = table(mmat_mex( kinf/2,x,signal)).Var1(1,1);
minfneg11 = table(mmat_mex(-kinf/2,x,signal)).Var1(1,1);
apos = 1/minfpos11;
aneg = 1/minfneg11;
N = round(angle(apos)-angle(aneg));
%[Kn,Cn]    = disdst(kz,x,signal,lambda,N);

%% Calculate analytical transfer matrix
M11t = zeros(size(M11));
M12t = zeros(size(M12));
bt   = zeros(size(M11));
i = 1;
for k = kz
    kap = k/2;
    xi = sqrt(lambda*eta0+kap^2);
    delta = 0.5*((xi/kap)+(kap/xi));
    gamma = 0.5*((xi/kap)-(kap/xi));
    M11_theory = (cos(2*xi*a) - 1i*delta*sin(2*xi*a))*exp(2*1i*kap*a);
    M12_theory = -1i*gamma*sin(2*xi*a);
    M11t(i) = (M11_theory);
    M12t(i) = (M12_theory);
    M21_theory = 1i*gamma*sin(2*xi*a);
    M22_theory = (cos(2*xi*a)+1i*delta*sin(2*xi*a))*exp(-2*1i*kap*a);
    M_theory = [M11_theory M12_theory; M21_theory M22_theory];
    b_theory = -M12_theory/M11_theory;
    bt(i) = b_theory;
    i = i+1;
end

figure; 
subplot(121);
plot(kz,M11);hold on; plot(kz, real(M11t)); legend('M11','M11t');
subplot(122);
plot(kz,M12);hold on; plot(kz, real(M12t)); legend('M12','M12t');

%% Calculate theoretical spectra components
kz_theory = linspace(0,0.2,5000)/(2*pi);
% 0.1 amp, positive wave -100<= x <= 100
fhatmag_theory = (0.1./(2*pi*kz_theory)).*abs(exp(-2*pi*1i*100*kz_theory) - exp(2*pi*1i*100*kz_theory)); 

%% Plot everything
figure;
subplot(151);
plot(x, signal);
xlabel('$x$','interpreter', 'latex');
title('Signal in Space Domain');

subplot(152);
semilogy(kz,abs(Fc)/(lambda*L),'.-');
xlim([0 0.2]);
%ylim([1e-12 1e1]);
title('Ordinary Fourier Transform');

subplot(153);
semilogy(kz,kz.*abs(b)/(lambda*L),'.-');
xlim([0 0.2]);
%ylim([1e-12 1e1]);
title('Direct Scattering Transform');

subplot(154);
semilogy(kz_theory*2*pi,fhatmag_theory,'.-');
xlim([0 0.2]);
%ylim([1e-12 1e1]);
title('Analytical Theory DFT');

subplot(155);
semilogy(kz,kz.*abs(bt)/(lambda*L),'.-');
xlim([0 0.2]);
%ylim([1e-12 1e1]);
title('Analytical Theory DFT');