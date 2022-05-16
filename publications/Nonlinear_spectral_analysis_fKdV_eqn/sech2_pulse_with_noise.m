% sech^2(x) pulse analytical
% Use data from 

clear all; close all; clc;
addpath('../../lib/')

L = 100; 
lambda = 1;
num_points = 150;
xz = linspace(-L/8, 7*L/8, num_points);
tend = 5;
num_points_t = 200;
tz = linspace(0, tend, num_points_t);
[XX,TT] = meshgrid(xz,tz);
c = 2;
x0 = 0;
full_solution = (c/2)*sech(sqrt(c)/2.*(XX-c*TT-x0)).^2;
%full_solution = sin(30*2*XX*pi/L); - for debugging FFT
lambdaD = 3.6e-2;
vIA = 5.5e3;
noise = 0.3*cos(2*pi*5e3*lambdaD/vIA.*TT) + 0.01*normrnd(0.0,0.1,num_points_t,num_points);
full_solution = full_solution + noise;

Emax = 5.0;
Emin = -2;
numE = 3000;

%% Perform calculation

tic
mu  = {};
fmu = {};
Ej  = {};
fEj = {};
N   = {};
Ez  = {};
M11 = {};
M12 = {};
Estar = {};
fEstar = {};
traceM = {};
m = {};
a = {}; % Amplitude of hyperelliptic functions

parfor s = 1:length(tz)    
    signal = full_solution(s,:);
    [mu_single, fmu_single, Ej_single, fEj_single, Estar_single, fEstar_single, N_single, Ez_single, M11_single,M12_single,traceM_single, m_single, a_single] = ist_spectra(Emin, Emax, numE, xz, signal*lambda);
    mu{s} = mu_single;
    fmu{s} = fmu_single;
    Ej{s} = Ej_single;
    fEj{s} = fEj_single;
    N{s} = N_single;
    Ez{s} = Ez_single;
    Estar{s} = Estar_single;        
    fEstar{s} = fEstar_single;
    M11{s} = M11_single;
    M12{s} = M12_single;
    traceM{s} = traceM_single;
    m{s} = m_single;
    a{s} = a_single;
    fprintf('Iteration %i / %i\n', s, length(tz));
end
toc

%% Post process IST stuff
num_modes = zeros(length(tz),1);
for s = 1:length(tz)
    num_modes(s) = N{s}-1;
end

num_solitons = zeros(length(tz),1);
k_ist = zeros(length(tz), max(num_modes));
a_ist = zeros(length(tz), max(num_modes));
m_ist = zeros(length(tz), max(num_modes));
for s = 1:length(tz)
    num_solitons(s) = length(find(m{s}>=0.999));
    k_ist(s, 1:num_modes(s)) = 2*pi*(1:num_modes(s))/L;
    a_ist(s, 1:num_modes(s)) = a{s};
    m_ist(s, 1:num_modes(s)) = m{s};
end

%% Compute FFT in space
Fs = num_points/L;            % Sampling frequency                    
T = 1/Fs;          % Sampling period     
Y = fft(full_solution,[],2); % FFT along columns (actually dim=2)
P2 = abs(Y/num_points);
P1 = P2(:,1:num_points/2+1);
P1(:,2:end-1) = 2*P1(:,2:end-1);
k = Fs*(0:(num_points/2))/num_points;

%% Ploting

figure;imagesc(tz,xz,full_solution');ylim([-5 15]);
set(gcf,'units','inches','position',[0,0,6,4]);
set(gca, 'FontName', 'Helvetica');
xlabel('t');
ylabel('x');
colorbar();
colormap(flipud(hot))
grid off;
box on;
ax = gca;
ax.LineWidth = 3;
saveas(gcf,'figures/homo_sech2_signal_with_noise.eps','epsc')
saveas(gcf,'figures/homo_sech2_signal_with_noise.png')

figure;imagesc(k,tz,P1);xlim([0 0.5]);
colorbar();
xlabel('k');
ylabel('t');
colormap("default")
set(gcf,'units','inches','position',[5,0,6,4]);
set(gca, 'FontName', 'Helvetica');
saveas(gcf,'figures/homo_sech2_fft_with_noise.eps','epsc')
saveas(gcf,'figures/homo_sech2_fft_with_noise.png')

figure;
plot(k_ist(1,:), a_ist(1,:),'bo');
ylabel('Amplitude of Mode')
xlabel('k')
hold on;
plot(k_ist(100,:), a_ist(100,:),'r*');
plot(k_ist(200,:), a_ist(200,:),'g.');
yyaxis right;
plot(k_ist(1,:), m_ist(1,:),'^');
plot(k_ist(100,:), m_ist(100,:),'>');
plot(k_ist(200,:), m_ist(200,:),'<');
legend(sprintf('t = %f', tz(1)),...
       sprintf('t = %f', tz(100)),...
       sprintf('t = %f', tz(200)),...
       sprintf('t = %f', tz(1)),...
       sprintf('t = %f', tz(100)),...
       sprintf('t = %f', tz(200)))
ylabel('Modulus');
xlim([0 0.5]);
set(gcf,'units','inches','position',[10,0,6,4]);
set(gca, 'FontName', 'Helvetica');
saveas(gcf,'figures/homo_sech2_ist_with_noise.eps','epsc')
saveas(gcf,'figures/homo_sech2_ist_with_noise.png')

%% Plot IST internal data
figure;
set(gcf,'units','inches','position',[0,4,4,8]);
set(gca, 'FontName', 'Cambria');
subplot(311);
xlim_n = [-0.75 0.8];
ylim_n = [-0.9 0.9];
idx = 2;
plot(Ez{idx}, M11{idx});
hold on;
title(sprintf('N = %i', N{idx}));
plot(Estar{idx}, fEstar{idx}, 'ro');symlog('y');
xlabel('E');
ylabel('$M_{11}$', 'interpreter', 'latex');
xlim(xlim_n)
ylim(ylim_n)
subplot(312);
plot(Ez{idx}, M12{idx}); 
hold on;
plot(mu{idx}, fmu{idx}, 'ro'); symlog('y');
xlabel('E');
ylabel('$M_{12}$', 'interpreter', 'latex');
xlim(xlim_n)
ylim(ylim_n)

subplot(313);
plot(Ez{idx}, traceM{idx});
hold on;
xlim(xlim_n)
ylim(ylim_n)

plot(Ej{idx}, fEj{idx}, 'ro');symlog('y');
xlabel('E');
ylabel('$\textit{tr} \textbf{M}/2$', 'interpreter', 'latex');
xlim(xlim_n)
ylim(ylim_n)
saveas(gcf,'figures/homo_sech2_ist_process_with_noise.eps','epsc')
saveas(gcf,'figures/homo_sech2_ist_process_with_noise.png')
