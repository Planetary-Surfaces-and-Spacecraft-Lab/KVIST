clear all; close all; clc;
% http://www1.spms.ntu.edu.sg/~ydchong/teaching/PH4401_Appendix_B_transfer_matrix.pdf

L = 1; 
lambda = 1e-2; 
num_points = 200;
xz = linspace(0, L, num_points);
delx = xz(2)-xz(1);
a = 000;
b = 200;
signal(xz<L/2) = a;
signal(xz>=L/2) = b;

Emax = 20.0;
Emin = -10.0;
Ez = linspace(Emin, Emax, 500);
kz = sqrt(Ez);
M11 = ones(length(kz), 1);
M11theory = ones(length(kz), 1);
M12 = zeros(length(kz), 1);
M12theory = zeros(length(kz), 1);
M21 = zeros(length(kz), 1);
M21theory = zeros(length(kz), 1);
M22 = ones(length(kz), 1);
M22theory = ones(length(kz), 1);
S = zeros(length(kz), 1);
traceM = ones(length(kz),1);
traceMtheory = ones(length(kz),1);
M11prev = 1;

i = 2;
for k = kz(2:end)
    [M,~,delx] = mmat_mex(Ez(i), xz, lambda*signal);
    
    %if(abs(det(M)-1.0) >= 3e-1 )
    %    error('oof, det(M) = %g, delx = %g, E = %g\n',det(M),delx,Ez(i));
    %end
    M11(i) = M(1,1);
    M12(i) = M(1,2);
    M21(i) = M(2,1);
    M22(i) = M(2,2);
    traceM(i) = 0.5*(M(1,1)+M(2,2));
    S(i) = S(i-1) + abs(sign(M11(i))-sign(M11prev))/2;
    M11prev = M11(i);
    
    % Calculate theoretical value of M from notes
    kapleft  = sqrt(lambda*a+Ez(i));
    kapright = sqrt(lambda*b+Ez(i));
    Mleft = [exp(1i*kapleft*L/2) 0; 0 exp(-1i*kapleft*L/2)];
    Mright = [exp(1i*kapright*L/2) 0; 0 exp(-1i*kapright*L/2)];
    Mstep  = 0.5*[1+(kapleft/kapright) 1-(kapleft/kapright); ...
                  1-(kapleft/kapright) 1+(kapleft/kapright)];
    Mtheory = Mright*Mstep*Mleft;
    Qleft = [1 1i*kapleft; 1 -1i*kapleft]; % Coordinate transformation matrix
                           % to change to identity basis 
                           % to compare apples to apples
    Qright = [1 1i*kapright; 1 -1i*kapright];
    Mtheory = inv(Qleft)*Mtheory*Qright;Mtheory=Mtheory';
    M11theory(i) = Mtheory(1,1);
    M21theory(i) = Mtheory(2,1);
    M22theory(i) = Mtheory(2,2);
    M12theory(i) = Mtheory(1,2);
    traceMtheory(i) = 0.5*trace(Mtheory);
    i = i+1;
    
end

figure;
subplot(231);
plot(Ez, M11,'.-');
xlim([Emin Emax]);
hold on;
plot(Ez, real(M11theory), '--');symlog('y');
xlabel('E');
yyaxis right;
plot(Ez, abs(M11-real(M11theory)));
ylabel('$M_{11}^{(n)} - M_{11}^{(t)}$','interpreter','latex');
legend('M_{11}', 'M_{11t}','error');

subplot(232);
plot(Ez, M12,'.-');
xlabel('E');
hold on;
plot(Ez, real(M12theory),'--');symlog('y');
xlim([Emin Emax]);
xlabel('E');
yyaxis right;
plot(Ez, abs(M12-real(M12theory)));
ylabel('$M_{12}^{(n)} - M_{12}^{(t)}$','interpreter','latex');
legend('M_{12}', 'M_{12t}','error');
xlim([Emin Emax]);

subplot(234);
plot(Ez, M21,'.-');
hold on;
plot(Ez, real(M21theory),'--');symlog('y');
xlabel('E');
yyaxis right;
plot(Ez, abs(M21-real(M21theory)));
ylabel('$M_{21}^{(n)} - M_{21}^{(t)}$','interpreter','latex');
legend('M_{21}', 'M_{21t}','error');
xlim([Emin Emax]);



subplot(235);
plot(Ez, M22,'.-');
hold on;
plot(Ez, real(M22theory),'--');symlog('y');
xlabel('E');
xlim([Emin Emax]);
yyaxis right;
plot(Ez, abs(M22-real(M22theory)));
ylabel('$M_{22}^{(n)} - M_{22}^{(t)}$','interpreter','latex');
legend('M_{22}', 'M_{22t}','error');

subplot(233);
plot(Ez, traceM,'.-');
hold on;
plot(Ez, real(traceMtheory),'--');symlog('y');
xlabel('E');
legend('Numerical', 'Theory');
xlim([Emin Emax]);
ylabel('$\frac{1}{2}tr M$','interpreter','latex');
subplot(236);
semilogy(Ez, abs(traceM-traceMtheory),'.-');
xlabel('E');
xlim([Emin Emax]);
ylabel('$\Delta \frac{1}{2}tr M$','interpreter','latex');