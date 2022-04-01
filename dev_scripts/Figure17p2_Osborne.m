% Benchmark against Osborne Chapt 17.6, pg 437

clear all; close all; clc;

L = 300; % Zab and Kru says this should be 2
lambda = 0.012; % Zab and Kru should be like 344.353
num_points = 300;
xz = linspace(0, L, num_points);
delx = xz(2)-xz(1);

% Initial conditions from 10.1103/PhysRevLett.15.240
signal = cos(pi*xz); % One full wave on [0,L], eqn 5 Z&K

fh = four(signal);
nyquistk = 2*pi/delx;
Emax = 0.06;
M11 = zeros(length(xz), 1);
M11(1) = 1;
S = zeros(length(xz),1);

i = 2;
for x = xz(i:end)
    [M,~] = mmat_mex(Emax, xz(1:i), lambda*signal(1:i));
    M11prev = M11(i-1);
    M11(i) = M(1,1);   
    S(i) = S(i-1) + abs(sign(M11(i))-sign(M11prev))/2;
    i = i+1;
end

figure;
plot(xz, M11);
hold on;
plot(xz, S);

% Matches figure 17.2 sort of. 
% Figure 17.2 seems to vary the frequency of the wave
% This seems to not do that

