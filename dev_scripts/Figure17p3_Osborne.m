% Benchmark against Osborne Chapt 17.3, pg 437
% Weird thing I just noticed (August 3, 2021):
%  Osborne uses units of cm^-1 in his book
%  they really ought to be cm^-2

clear all; close all; clc;

L = 300; % Zab and Kru says this should be 2
%lambda = 0.012; % Zab and Kru should be like 344.353
lambda = 0.012;
num_points = 300;
xz = linspace(0, L, num_points);
delx = xz(2)-xz(1);

% Initial conditions from 10.1103/PhysRevLett.15.240
signal = cos(pi*xz); % One full wave on [0,L], eqn 5 Z&K

Emax = 0.06;
Emin = -0.02;
Ez = linspace(Emin, Emax, 1000);
kz = sqrt(Ez);
M11 = ones(length(kz), 1);
S = zeros(length(kz), 1);
traceM = ones(length(kz),1);
M11prev = 1;

i = 2;
for k = kz(2:end)
    [M,~] = mmat_mex(Ez(i), xz, lambda*signal);
    M11(i) = M(1,1);
    traceM(i) = 0.5*(M(1,1)+M(2,2));
    S(i) = S(i-1) + abs(sign(M11(i))-sign(M11prev))/2;
    M11prev = M11(i);
    i = i+1;
end

figure;
subplot(121);
plot(Ez, real(traceM),'.-');
ylim([-15, 20]);
xlim([Emin Emax]);
hold on;
yyaxis right;
plot(Ez, S);
ylim([0 28]);
subplot(122);
plot(Ez, real(traceM),'.');
ylim([-5, 7]);
xlim([Emin Emax]);
hold on;
yyaxis right;
plot(Ez, S);
xlim([0 0.024]);
ylim([8 17]);

figure;
subplot(121);
plot(Ezcondst, real(traceMcondst)/2,'.-');
ylim([-15, 20]);
xlim([Emin Emax]);
hold on;
yyaxis right;
plot(Ez, S);
ylim([0 28]);
subplot(122);
plot(Ez, real(traceMcondst)/2,'.');
ylim([-5, 7]);
xlim([Emin Emax]);
hold on;
yyaxis right;
plot(Ez, S);
xlim([0 0.024]);
ylim([8 17]);