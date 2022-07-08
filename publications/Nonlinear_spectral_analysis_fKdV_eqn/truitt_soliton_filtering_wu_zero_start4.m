clear all; close all; clc;
addpath('../../lib/')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Wu forced zero start test case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = readmatrix('truitt_data/wu_zero_start/U_wu_zero_start.txt');
xz= readmatrix('truitt_data/wu_zero_start/xz_wu_zero_start.txt');
tz = readmatrix('truitt_data/wu_zero_start/tz_wu_zero_start.txt');


%% Clip data to interesting part
Emax = 2;
Emin = -1.4;
numE = 50000;
alpha = -3/2;
beta = -1/6;
lambda = alpha/(6*beta);
signal = real(U(800:1250,40))';
xz_original = xz;
xz = xz(800:1250);


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

E1 = Ej{1}(1);

%% Post process IST stuff
L = xz(end)-xz(1);
num_modes = zeros(length(tz),1);
for s = 1:length(tz)
    num_modes(s) = N{s}-1;
end

num_solitons = zeros(length(tz),1);
k_ist = zeros(length(tz), max(num_modes));
a_ist = zeros(length(tz), max(num_modes));
for s = 1:length(tz)
    num_solitons(s) = length(find(m{s}>=0.999));
    k_ist(s, 1:num_modes(s)) = 2*pi*(1:num_modes(s))/L;
    a_ist(s, 1:num_modes(s)) = a{s};
end

%% Ploting
% .mat file for above analysis: truitt_data/wu_zero_start/ist_analysis4.mat

figure;
imagesc(tz,flip(xz_original),flip(real(U)));
set(gcf,'units','inches','position',[0,0,9,4]);
set(gca, 'FontName', 'Arial')
colorbar();
colormap('jet');
xlabel('t');
ylabel('x');
yline(xz(1),'r','linewidth',2);
yline(xz(end),'r','linewidth',2);
xline(tz(20),'c','linewidth',2);
xline(tz(55),'m','linewidth',2);
xline(tz(90),'g','linewidth',2);
xline(tz(40),'o','linewidth',2);
ylim([-50 50]);
saveas(gcf,'figures/wu_zero_sim_overview.eps','epsc')
saveas(gcf,'figures/wu_zero_sim_overview.png')

figure; 
set(gcf,'units','inches','position',[5,0,4,4]);
set(gca, 'FontName', 'Helvetica')
plot(xz, signal);
xlabel('x');
grid on;
saveas(gcf,'figures/wu_zero_start_signal4_no_ylabel.eps','epsc')
saveas(gcf,'figures/wu_zero_start_signal4_no_ylabel.png');
ylabel('u(x)');
saveas(gcf,'figures/wu_zero_start_signal4.eps','epsc')
saveas(gcf,'figures/wu_zero_start_signal4.png')


figure;
set(gcf,'units','inches','position',[10,0,10,4]);
set(gca, 'FontName', 'Arial');
for i = 1:6
    subplot(320+i);
    plot(xz,modemat(i,:));
    ylabel(sprintf("u_%i(x)", i));
    if i == 5 || i == 6
        xlabel('x');
    end
end
saveas(gcf,'figures/wu_zero_start_decomp4.eps','epsc')
saveas(gcf,'figures/wu_zero_start_decomp4.png')

%%
figure;
plot(k_ist(1,:), a_ist(1,:),'o');
ylabel('Amplitude of Mode')
xlabel('k')
hold on;
yyaxis right;
plot(k_ist(1,:), nine_digcnt(m{1}),'^');
hold on;
yline(3,'--','LineWidth',2, ...
      'Color',[0.8500 0.3250 0.0980]);
ylabel('\sigma (Modulus)');
xlim([0 2.0]);
legend('Amplitude', '\sigma (Modulus)');
grid minor;
set(gcf,'units','inches','position',[8,5,5,4]);
set(gca, 'FontName', 'Helvetica');
saveas(gcf,'figures/wu_zero_start_ist4.eps','epsc')
saveas(gcf,'figures/wu_zero_start_ist4.png')
