clear all; close all; clc;
addpath('../../lib/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Giak test case (single pinned soliton) 12288 cells in space axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% U = readmatrix('truitt_data/gaik_12288/gaik_12288_U.txt');
% xz= readmatrix('truitt_data/gaik_12288/gaik_12288_x.txt');
% tz = readmatrix('truitt_data/gaik_12288/gaik_12288_t.txt');
U = readmatrix('truitt_data/gaik_time_series_1000/gaik_1000_tz_U.txt');
xz= readmatrix('truitt_data/gaik_time_series_1000/gaik_1000_tz_x.txt');
tz = readmatrix('truitt_data/gaik_time_series_1000/gaik_1000_tz_t.txt');

Emax = 5;
Emin = -11;
numE = 700;
alpha = 6;
beta = 1;
lambda = alpha/(6*beta);
signal = real(U(6000:6300,1))';
xz_sim = xz;
xz = xz(6000:6300);


figure;
plot(xz,signal/max(signal));

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

for s = 1:length(xz)    
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

E1 = Ej{1}(1);

%% Ploting

figure;imagesc(tz,xz_sim,real(U));ylim([-5 130])
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
saveas(gcf,'figures/homo_three_soliton_signal.eps','epsc')
saveas(gcf,'figures/homo_three_soliton_signal.png')

%% More ploting
figure; 
set(gcf,'units','inches','position',[0,0,5,4]);
set(gca, 'FontName', 'Arial')
plot(xz, signal);
xlabel('x');
ylabel('u(x)');
grid on;
saveas(gcf,'figures/homo_three_soliton_initial_signal.eps')
saveas(gcf,'figures/homo_three_soliton_initial_signal.png')


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
saveas(gcf,'figures/homo_three_soliton_decomp.eps')
saveas(gcf,'figures/homo_three_soliton_decomp.png')


