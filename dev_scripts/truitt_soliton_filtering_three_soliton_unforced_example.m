clear all; close all; clc;

%% DONT CHANGE MY PARAMETERS - I WORK AS IS


% load 'truitt_2019_solition_solution.mat';
% % Interesting locations, 400, 1750
% signal = real(U(1750,:));
% xz = tvec;
% L = xz(end)-xz(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Giak test case (single pinned soliton) 12288 cells in space axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% U = readmatrix('truitt_data/gaik_12288/gaik_12288_U.txt');
% xz= readmatrix('truitt_data/gaik_12288/gaik_12288_x.txt');
% tz = readmatrix('truitt_data/gaik_12288/gaik_12288_t.txt');
U = readmatrix('truitt_data/gaik_time_series_1000/gaik_1000_tz_U.txt');
xz= readmatrix('truitt_data/gaik_time_series_1000/gaik_1000_tz_x.txt');
tz = readmatrix('truitt_data/gaik_time_series_1000/gaik_1000_tz_t.txt');


c0 = 1; % Moving reference frame speed

% Change to stationary reference frame
tform = affine2d([1 c0 0; 0 1 0; 0 0 1]);
[TZ, XZ] = meshgrid(tz,xz);
U_stationary = imwarp(U, tform);
XZ_stationary = imwarp(XZ, tform);
TZ_stationary = imwarp(TZ, tform);
dx = xz(2)-xz(1);
xz_stationary = [xz xz(end)+dx:dx:(length(XZ_stationary)-length(XZ))*dx+xz(end)]';

plot_comparison_of_sheared_images = 0;
if (plot_comparison_of_sheared_images)
   figure;
   subplot(121);
   image(xz,tz,100*abs(U'));
   set(gca,'YDir','normal');
   xlabel('x');
   ylabel('t');
   title('U moving');
   subplot(122);
   image(xz_stationary,tz,100*abs(U_stationary'));
   set(gca,'YDir','normal') ;
   title('U still');
   xlabel('x');
   ylabel('t');
   
   figure;
   subplot(121);
   image(xz,tz,abs(XZ)');
   set(gca,'YDir','normal');
   title('XZ moving');
   xlabel('x');
   ylabel('t');
   subplot(122);
   image(xz_stationary,tz,abs(XZ_stationary)');
   set(gca,'YDir','normal');
   title('XZ still');
   xlabel('x');
   ylabel('t');
   
   figure;
   subplot(121);
   image(xz,tz,100*abs(TZ)');
   set(gca,'YDir','normal');
   xlabel('x');
   ylabel('t');
   title('TZ moving');
   subplot(122);
   image(xz_stationary,tz,100*abs(TZ_stationary)');
   set(gca,'YDir','normal');
   title('TZ still');
   xlabel('x');
   ylabel('t');
end


% 
% %%
% % Shift moving reference frame back to still reference frame
% waitb=waitbar(0.0, 'Moving to Stationary Frame');
% for xi = 1:length(xz_stationary)
%     waitbar(xi/length(xz_stationary), waitb, 'Moving to Stationary Frame');
%     for ti = 1:length(tz_stationary)
%         x = xz_stationary(xi);
%         t = tz_stationary(ti);
%         xmoving = x-c0*t;
%         if ismember(xmoving,xz) 
%             movingidx = find(xz==xmoving);
%             U_stationary(xi, ti) = U(movingidx,ti);
%         else
%             U_stationary(xi, ti) = nan; % no data exists here
%         end
%     end
% end

%%
% Downsample in spatial domain
% U = downsample(U,12);
% xz = downsample(xz,12);


% WONT RECOGNIZE A SOLITON IN TIME IF SOLUTIONS ARE IN SPACE
% WIDTH OF SOLITON IN TIME DOES NOT EQUAL HEIGHT

L = xz(5900:6400);
%L = tz(850)-tz(750);

Emax = 10;
Emin = -15;
numE = 300;
alpha = 6;
beta = 1;
lambda = alpha/(6*beta);
signal = real(U(5900:6400,1))';
%signal = real(U_stationary(8000,500:600));
%xz=tz(500:600);
xz = xz(5900:6400);

% fitted signal
center = 1.707;
c = 4.0*10000;
signal_fit = (c/2)*sech(sqrt(c/2)*(xz-center)).^2;

figure;
plot(xz,signal/max(signal));
hold on;
plot(xz, signal_fit/c*2);
legend('Data', 'Fit'); 

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

% TODO: LOOK AT M12 PASS. WHY IS IT JUMPING TO RIGHT BOUNDARY????!!!!!
%       ADD IN DEBUGGING INFO FOR CUSTOM INTERVAL SEARCH
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
% subplot(623);
% plot(xz, modemat(2,:));
% xlim([xz(1) xz(end)]);
% xlabel('x');
% ylabel('Mode 2');
% subplot(624);
% plot(xz, modemat(3,:));
% xlim([xz(1) xz(end)]);
% xlabel('x');
% ylabel('Mode 3');
% subplot(625);
% plot(xz, modemat(4,:));
% xlim([xz(1) xz(end)]);
% xlabel('x');
% ylabel('Mode 4');
% subplot(626);
% plot(xz, modemat(5,:));
% xlim([xz(1) xz(end)]);
% xlabel('x');
% ylabel('Mode 5');
% % subplot(627);
% plot(xz, modemat(6,:));
% xlim([xz(1) xz(end)]);
% xlabel('x');
% ylabel('Mode 6');
% subplot(627);
% plot(xz, modemat(6,:));
% xlim([xz(1) xz(end)]);
% xlabel('x');
% ylabel('Mode 6');
% subplot(628);
% plot(xz, modemat(7,:));
% xlim([xz(1) xz(end)]);
% xlabel('x');
% ylabel('Mode 7');
% subplot(629);
% plot(xz, modemat(8,:));
% xlim([xz(1) xz(end)]);
% xlabel('x');
% ylabel('Mode 8');
% subplot(6,2,10);
% plot(xz, modemat(9,:));
% xlim([xz(1) xz(end)]);
% xlabel('x');
% ylabel('Mode 9');
% subplot(6,2,11);
% plot(xz, modemat(10,:));
% xlim([xz(1) xz(end)]);
% xlabel('x');
% ylabel('Mode 10');
% subplot(6,2,12);
% plot(xz, modemat(11,:));
% xlim([xz(1) xz(end)]);
% xlabel('x');
% ylabel('Mode 11');
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


