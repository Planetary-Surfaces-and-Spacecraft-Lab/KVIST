function plot_M(mu, fmu, Ej, fEj, Estar, fEstar, N, Ez, M11,M12,traceM, m, a)

figure;
subplot(311);
plot(Ez, M11);
hold on;
title(sprintf('N = %i', N));
plot(Estar, fEstar, 'ro');symlog('y');
legend('M_{11}', 'E^{*}');
xlabel('E');
ylabel('$M_{11}$', 'interpreter', 'latex');
subplot(312);
plot(Ez, M12); 
hold on;
plot(mu, fmu, 'ro'); symlog('y');
xlabel('E');
ylabel('$M_{12}$', 'interpreter', 'latex');
xline(Estar, 'm--');
legend('M_{12}', '\mu_j', 'E^{*}');

subplot(313);
plot(Ez, traceM);
hold on;
plot(Ej, fEj, 'ro');symlog('y');
legend('tr M /2', 'E_j');
xline(mu, 'm--');
xlabel('E');
ylabel('tr M /2');
legend('tr M /2', 'E_j', '\mu');
end