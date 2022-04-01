clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Debris forcing solitons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = readmatrix('truitt_data/DebrisD/U_DebrisD.txt');
xz= readmatrix('truitt_data/DebrisD/xz_DebrisD.txt');
tz = readmatrix('truitt_data/DebrisD/tz_DebrisD.txt');

Emax = 5;
Emin = -0.2;
numE = 3000;
alpha = 1.0051;
beta = 0.4925;
lambda = alpha/(6*beta);
signal = real(U(3500:4200,50))';
signal = signal+randn(1,length(signal))*max(signal); % 100 percent white noise
xz = xz(3500:4200);


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

figure;subplot(621);
plot(xz, E1*ones(size(xz)));
ylabel('$E_1$','interpreter', 'latex');
xlim([xz(1) xz(end)]);
xlabel('x');
subplot(622);
plot(xz, modemat(1,:));
xlabel('x');
xlim([xz(1) xz(end)]);
ylabel('Mode 1');
subplot(623);
plot(xz, modemat(2,:));
xlim([xz(1) xz(end)]);
xlabel('x');
ylabel('Mode 2');
subplot(624);
plot(xz, modemat(3,:));
xlim([xz(1) xz(end)]);
xlabel('x');
ylabel('Mode 3');
subplot(625);
plot(xz, modemat(4,:));
xlim([xz(1) xz(end)]);
xlabel('x');
ylabel('Mode 4');
subplot(626);
plot(xz, modemat(5,:));
xlim([xz(1) xz(end)]);
xlabel('x');
ylabel('Mode 5');
subplot(627);
plot(xz, modemat(6,:));
xlim([xz(1) xz(end)]);
xlabel('x');
ylabel('Mode 6');
subplot(627);
plot(xz, modemat(6,:));
xlim([xz(1) xz(end)]);
xlabel('x');
ylabel('Mode 6');
subplot(628);
plot(xz, modemat(7,:));
xlim([xz(1) xz(end)]);
xlabel('x');
ylabel('Mode 7');
subplot(629);
plot(xz, modemat(8,:));
xlim([xz(1) xz(end)]);
xlabel('x');
ylabel('Mode 8');
subplot(6,2,10);
plot(xz, modemat(9,:));
xlim([xz(1) xz(end)]);
xlabel('x');
ylabel('Mode 9');
subplot(6,2,11);
plot(xz, modemat(10,:));
xlim([xz(1) xz(end)]);
xlabel('x');
ylabel('Mode 10');
subplot(6,2,12);
plot(xz, modemat(11,:));
xlim([xz(1) xz(end)]);
xlabel('x');
ylabel('Mode 11');
% figure;
% plot(Ez,M12/10,'k--');symlog('y');%,'DisplayName','M_{12}');
% hold on;
% plot(Ez,real(traceM),'k');%,'DisplayName','tr M/2');
% %plot(Ez,M11,'c--');symlog('y');
% plot(Ej, fEj, 'bo');symlog('y');%,'DisplayName','E_j');
% plot(mu,fmu,'ro');symlog('y');%,'DisplayName','\mu');
% %plot(Estar,fEstar,'go');symlog('y');%,'DisplayName','E^*');symlog('y');
% legend('M_{12}','tr M /2','E_j','\mu');
% ylim([-3, 3]);
% xlim([Ez(1) Ez(end)]);

figure;
subplot(131);
plot(xz, signal);
title('Original Signal');
ylabel('u(x)');
subplot(132);
if (lambda ~= 0.0)
    plot(xz, ((-E1+sum(modemat,1))/lambda));
else
    plot(xz, ((-E1+sum(modemat,1))));
end
title('Reconstructed Signal');
xlabel('x');
subplot(133);
xlabel('x');
if (lambda ~= 0.0)
    semilogy(xz, abs((signal)-(-E1+sum(modemat,1)))/lambda)
else
    semilogy(xz, abs((signal)-(-E1+sum(modemat,1))))
end
title('Total Error');

% denoising:
figure; 
noise = signal - real(U(3500:4200,50))';
subplot(241);
plot(xz,real(U(3500:4200,50))');
title('True signal')
ylim([-0.4 0.6]);
xlim([min(xz) max(xz)]);
subplot(242);
plot(xz, ((-E1+sum(modemat(1,:),1))/lambda))
title('Mode 1')
ylim([-0.4 0.6]);
xlim([min(xz) max(xz)]);
subplot(243);
plot(xz, ((-E1+sum(modemat(2,:),1))/lambda))
title('Mode 2')
ylim([-0.4 0.6]);
xlim([min(xz) max(xz)]);
subplot(244);
plot(xz, ((-E1+sum(modemat(8,:),1))/lambda))
title('Mode 8')
ylim([-0.4 0.6]);
xlim([min(xz) max(xz)]);
subplot(245);
plot(xz,noise);
title('Noise')
ylim([-0.4 0.6]);
xlim([min(xz) max(xz)]);
subplot(246);
plot(xz,((-E1+sum(modemat(2:4,:),1))/lambda));
title('Modes 2:4')
ylim([-0.4 0.6]);
xlim([min(xz) max(xz)]);
subplot(247);
plot(xz,((-E1+sum(modemat(5:end,:),1))/lambda));
title('Modes 5:end')
ylim([-0.4 0.6]);
xlim([min(xz) max(xz)]);
subplot(248);
plot(xz,((-E1+sum(modemat(1:5,:),1))/lambda));
title('Modes 1:5')
ylim([-0.4 0.6]);
xlim([min(xz) max(xz)]);

figure; 
subplot(311)
plot(xz,real(U(3500:4200,50))');
ylabel('U')
ylim([-0.3 0.5]);
xlim([min(xz) max(xz)]);
text(0.9*max(xz),0.35,'(a)','FontSize',16)
title('Simulation signal');
subplot(312)
plot(xz,((-E1+sum(modemat(1:6,:),1))/lambda));
ylabel('U')
ylim([-0.3 0.5]);
text(0.9*max(xz),0.35,'(b)','FontSize',16)
xlim([min(xz) max(xz)]);
title('Filtered signal');
subplot(313)
plot(xz,signal)
xlim([min(xz) max(xz)]);
text(0.9*max(xz),1.4,'(c)','FontSize',16)
xlabel('X'); ylabel('U')
title('Measured signal');


