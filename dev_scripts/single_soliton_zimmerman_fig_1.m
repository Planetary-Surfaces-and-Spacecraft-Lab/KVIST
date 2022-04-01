% https://npg.copernicus.org/articles/6/11/1999/npg-6-11-1999.pdf
% https://doi.org/10.5194/npg-6-11-1999
% Figure 1

clear all; close all; clc;

% August 17, this seems to work better with 
% cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 2, 2, 1.0, T, 2, M, 2, 0.0, Mtmp, 2);
% Although the negative stuff still looks whack

% August 24, det T = 1, which is good. I though cosh(x) and sinh(x) were
% getting outside the range of numerical accuracy, but that doesn't seem to
% be the case (it still can happen for very negative kappa*delx). 
% Works for positive kap*delx which seems to indicate the MATLAB bindings are
% working. 

L = 10; 
lambda = 1;
c= 2;
num_points = 500;
xz = linspace(0, L, num_points);
signal = (c/2)*sech(sqrt(c)/2*(xz-(L/2))).^2;

Emax = 3.0;
Emin = -1.5;
Ez = linspace(Emin, Emax, 1000);
kz = sqrt(Ez);
M11 = ones(length(kz), 1);
M12 = zeros(length(kz), 1);
S = zeros(length(kz), 1);
traceM = ones(length(kz),1);
M11prev = 1;

i = 2;
for k = kz(2:end)
    [M,~] = mmat_mex(Ez(i), xz, lambda*signal);
    if(abs(det(M)-1.0) >= 3e-2 )
        error('oof');
    end
    M11(i) = M(1,1);
    M12(i) = M(1,2);
    traceM(i) = 0.5*(M(1,1)+M(2,2));
    S(i) = S(i-1) + abs(sign(M11(i))-sign(M11prev))/2;
    M11prev = M11(i);
    i = i+1;
    
end

figure;
subplot(311);
plot(Ez, traceM,'.-');
hold on;
plot(Ez, ones(size(Ez)), 'k--', 'LineWidth',0.5);symlog('y');
plot(Ez, -ones(size(Ez)), 'k--', 'LineWidth',0.5);symlog('y');
ylim([-2, 5]);
xlim([Emin Emax]);
title('tr M/2');
xlabel('E');
hold on;
yyaxis right;
plot(Ez, S);
subplot(312);
plot(Ez, M12,'.-');
hold on;
plot(Ez, ones(size(Ez)), 'k--', 'LineWidth',0.5);symlog('y');
plot(Ez, -ones(size(Ez)), 'k--', 'LineWidth',0.5);symlog('y');
xlabel('E');
ylim([-2, 5]);
title('M_{12}');
xlim([Emin Emax]);
hold on;
yyaxis right;
plot(Ez, S);
ylim([8 17]);
subplot(313);
plot(xz, signal);
xlabel('x');
ylabel('u(x)');
title('signal');