clear all; clc; close all; 
addpath('../../lib/')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Debris forcing solitons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = readmatrix('truitt_data/DebrisD/U_DebrisD.txt');
xz= readmatrix('truitt_data/DebrisD/xz_DebrisD.txt');
tz = readmatrix('truitt_data/DebrisD/tz_DebrisD.txt');

% Conditions in paper
Emax = 3;
Emin = -1.5;
numE = 3000;

%tz_PIST_idx = [45 65 90];
tz_PIST_idx = [65]; % was 66 when I was talking to Christine
tz_PIST = tz(tz_PIST_idx);
alpha = 1.051;
beta = 0.4925;
lambda = alpha/(6*beta);
%xclip_min_idx = 4055; % both solitons (moving)
xclip_min_idx = 4235; % second soliton clipped out (moving)
%xclip_min_idx = 3955; % both solitons (debris not moving)
%xclip_min_idx = 4099; % second soliton clipped out (debris not moving)
%xclip_min_idx = 4005;
xclip_max_idx = 4505;
signal = real(U(xclip_min_idx:xclip_max_idx,:))';
% signal = [signal; flip(signal)];
% signal = signal(floor(length(signal)/4):floor(3*length(signal)/4)-1,:)';
xz_original = xz;
xz = xz(xclip_min_idx:xclip_max_idx);

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

for s = 1:length(tz_PIST_idx)
    [mu_single, fmu_single, Ej_single, fEj_single, Estar_single, ...
     fEstar_single, N_single, Ez_single, M11_single,M12_single,...
     traceM_single, m_single, a_single] = ...
              ist_spectra(Emin, Emax, numE, xz, signal(tz_PIST_idx(s),:));
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
    fprintf('Iteration %i / %i\n', s, length(tz_PIST_idx));
end
toc


%% Plotting

% Plot entire time series
figure;
imagesc(tz,xz,signal');
set(gcf,'units','inches','position',[0,0,8,4]);
set(gca, 'FontName', 'Helvetica')
grid off;
box on;
ax = gca;
ax.LineWidth = 2;
colorbar();
colormap(flipud(hot))
xlabel('t');
ylabel('x');
%yline(xz(1),'r','linewidth',2);
%yline(xz(end),'r','linewidth',2);
lwl = 3;
for i = tz_PIST_idx
    xline(tz(i),'b','linewidth',lwl);
end
saveas(gcf,'figures/debris_forcing_sim_clipped_for_second_soliton.eps','epsc')
saveas(gcf,'figures/debris_forcing_sim_clipped_for_second_soliton.png')

for i = 1:length(tz_PIST_idx)
    figure;
    num_modes = length(a{i});
    L = max(xz) - min(xz);
    plot(2*pi*(1:num_modes)/L, a{i},'o');
    ylabel('Amplitude of Mode')
    xlabel('k')
    hold on;
    yyaxis right;
    plot(2*pi*(1:num_modes)/L, nine_digcnt(m{i}),'^');
    hold on;
    plot(2*pi*(0:num_modes)/L,3*ones(length(2*pi*(0:num_modes)/L)),'--','LineWidth',2);
    ylabel('\sigma (Modulus)');
    xlim([0 2.0]);
    legend('Amplitude', '\sigma (Modulus)');
    grid minor;
    set(gcf,'units','inches','position',[8,5,5,4]);
    set(gca, 'FontName', 'Helvetica');
    saveas(gcf,sprintf('figures/debris_forcing_ist%i_clipped_for_second_soliton.eps',i),'epsc')
    saveas(gcf,sprintf('figures/debris_forcing_ist%i_clipped_for_second_soliton.png',i))
end





