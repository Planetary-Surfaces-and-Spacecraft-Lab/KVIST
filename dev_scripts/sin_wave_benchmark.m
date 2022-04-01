% Benchmark against Osborne Chapt 17.6, pg 437
% Zabusky and Kruskal problem

clear all; %close all; clc;
clc;

L = 2*pi; % Zab and Kru says this should be 2
lambda = 10;%0.022; % Zab and Kru should be like 344.353
num_points = 1000;
xz = linspace(0, L, num_points);
delx = xz(2)-xz(1);

% Initial conditions from 10.1103/PhysRevLett.15.240
signal = 1*sin(xz); % One full wave on [0,L], eqn 5 Z&K

% 09/16/21 - MATCHES ZIMMERMAN & HAARLEMMER FIG 2

Emax = 40;
Emin = -20;
Ez = linspace(Emin, Emax, 2000);
kz = sqrt(Ez);
M11 = ones(length(kz), 1);
M12 = zeros(length(kz), 1);
S = zeros(length(kz), 1);
traceM = ones(length(kz),1);
M11prev = 1;

i = 2;
for k = kz(2:end)
    [M,~] = mmat_mex(Ez(i), xz, lambda*signal);
%     if(abs(det(M)-1.0) >= 3e-2 )
%         error('oof');
%     end
    M11(i) = M(1,1);
    M12(i) = M(1,2);
    traceM(i) = 0.5*(M(1,1)+M(2,2));
    S(i) = S(i-1) + abs(sign(M11(i))-sign(M11prev))/2;
    M11prev = M11(i);
    i = i+1;
    
end

%kinf = 1e18;
%minfpos11 = table(mmat( kinf/2,x,signal,lambda)).Var1( 1,1);
%minfneg11 = table(mmat(-kinf/2,x,signal,lambda)).Var1(1,1);
%apos = 1/minfpos11;
%aneg = 1/minfneg11;
%N = round(angle(apos)-angle(aneg));
%[Kn,Cn]    = disdst(kz,x,signal,lambda,N);

figure('Position',[600 300 1600,800]);
subplot(221);
plot(xz, signal);
xlabel('$x$','interpreter', 'latex');
title('Signal in Space Domain');
% 
% subplot(132);
% four_amp = abs(fh);
% semilogy(kz,abs(four_amp(1:501))/(lambda*L),'o-');
% xlim([0 2*5*2*pi/L]);
% %xlim([kz(1) kz(round(length(kz)/4))]);
% ylim([1e-12 1e1]);
% title('Ordinary Fourier Transform');
% 
% subplot(133);
% % semilogy(kz,kz.*abs(b)/(lambda*L),'o-');
% semilogy(kz,abs(Fc)/(lambda*L),'o-');
% xlim([0 2*5*2*pi/L]);
% %xlim([kz(1) kz(round(length(kz)/4))]);
% ylim([1e-12 1e1]);
% title('Direct Scattering Transform Comp $F_c$');

subplot(222);
plot(Ez(2:end),traceM(2:end));symlog('y');
xlabel('E');
ylabel('tr M / 2');

subplot(223);
plot(Ez(2:end),M11(2:end));symlog('y');
xlabel('E');
ylabel('M_{11}');

subplot(224);
plot(Ez(2:end),M12(2:end));symlog('y');
xlabel('E');
ylabel('M_{12}');