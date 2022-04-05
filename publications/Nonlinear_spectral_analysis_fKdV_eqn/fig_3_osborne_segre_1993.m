% "A numerical inverse scattering transform for the periodic KdV equation"
% Figure 3

clear all; close all; clc;
addpath('../../lib/')

L = 10; 
lambda = 1;
c= 2;
num_points = 100;
xz = linspace(0, L, num_points);
signal = 2.0*gaussmf(xz,[1.0,5]);

Emax = 5.0;
Emin = -2.0;
numE = 200; % 1000

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

parfor s = 1:length(xz)    
    [mu_single, fmu_single, Ej_single, fEj_single, Estar_single, fEstar_single, N_single, Ez_single, M11_single,M12_single,traceM_single, m_single, a_single] = ist_spectra(Emin, Emax, numE, xz, circshift(signal,-(s-1))*lambda);
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
    fprintf('Iteration %i / %i\n', s, length(xz));
end
toc

mumat = zeros(N{1}-1,length(xz));
for j = 1:length(xz)
    if length(mumat(:,j)) == length(mu{j})
        mumat(:,j) = mu{j};
    else
        fprintf('j, err = %i\n', j);
    %    mumat(:,j) = padarray(mu{j}, N{1}-1-length(mu{j}),'post');
    end
end

modemat = zeros(N{1}-1,length(xz));
for j = 1:(N{1}-1)
    modemat(j,:) = (2*mumat(j,:) - Ej{1}(2*j) - Ej{1}(2*j+1)); 
end

%% Ploting

figure; 
set(gcf,'units','inches','position',[0,0,5,4]);
set(gca, 'FontName', 'Arial')
plot(xz, signal);
xlabel('x');
ylabel('u(x)');
grid on;
saveas(gcf,'figures/gaussian_pulse_signal.eps')
saveas(gcf,'figures/gaussian_pulse_signal.png')


figure;
set(gcf,'units','inches','position',[0,0,10,4]);
set(gca, 'FontName', 'Arial');
for i = 1:6
    subplot(320+i);
    plot(xz,modemat(i,:));
    ylabel(sprintf("u_%i(x)", i));
    if i == 5 || i == 6
        xlabel('x');
    end
end
saveas(gcf,'figures/gaussian_pulse_decomp.eps')
saveas(gcf,'figures/gaussian_pulse_decomp.png')


