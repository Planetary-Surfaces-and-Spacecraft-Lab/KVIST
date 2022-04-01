% "A numerical inverse scattering transform for the periodic KdV equation"
% Figure 3

clear all; close all; clc;

L = 10; 
lambda = 1;
c= 2;
num_points = 100;
xz = linspace(0, L, num_points);
signal = 2.0*gaussmf(xz,[1.0,5]);

Emax = 5.0;
Emin = -2.0;
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
subplot(221);
plot(Ez, traceM,'.-');symlog('y');
xlabel('E');
hold on;
plot(Ez, ones(size(Ez)), 'k--', 'LineWidth',0.5);symlog('y');
plot(Ez, -ones(size(Ez)), 'k--', 'LineWidth',0.5);symlog('y');
xlim([Emin Emax]);
ylim([-3 5]);
title('tr M');
hold on;
yyaxis right;
plot(Ez, S);
xlabel('E');
subplot(223);
plot(Ez, M12,'.-');
hold on;
plot(Ez, ones(size(Ez)), 'k--', 'LineWidth',0.5);
plot(Ez, -ones(size(Ez)), 'k--', 'LineWidth',0.5);symlog('y');
xlim([Emin Emax]);
ylim([-3 5]);
xlabel('E');
title('M_{12}');
hold on;
yyaxis right;
plot(Ez, S);
ylim([8 17]);
subplot(222);
plot(xz, signal);
xlabel('x');
ylabel('u(x)');
title('signal');