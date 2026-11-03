clear;
clc;
close all;

source('parte1.m');

% Datos del problema
L = 1;
alpha = 0.01;

u0 = @(x) 80*sin(pi*x/L);

T_exacta = @(x,t) 20 + 80*exp(-alpha*(pi/L)^2*t).*sin(pi*x/L);

% Experimento base
M = 50;
dt = 0.01;
tf = 10;

dx = L/M;
mu = alpha*dt/dx^2;

tic;

[x_exp,t_exp,U_exp] = euler_explicito(L,alpha,u0,M,dt,tf);

T_exp = U_exp + 20;

tiempo_exp = toc;

% Verificar estabilidad
fprintf('\nEuler explicito\n');
fprintf('mu = %.4f\n',mu);

if mu <= 0.5
    fprintf('Condicion de estabilidad: cumple\n');
else
    fprintf('Condicion de estabilidad: no cumple\n');
end

% Error al tiempo final
T_ex_final = T_exacta(x_exp,tf);

E_f_exp = max(abs(T_exp(end,:) - T_ex_final));

% Error en toda la malla
E_malla_exp = 0;

for n = 1:length(t_exp)

    T_ex_n = T_exacta(x_exp,t_exp(n));

    error_n = max(abs(T_exp(n,:) - T_ex_n));

    E_malla_exp = max(E_malla_exp,error_n);

end

% Resultados
fprintf('\n');
fprintf('%-18s %-14s %-14s %-12s %-12s %-12s\n','Metodo','E_f','E_malla','Tiempo(s)','Residuo','Iteraciones');
fprintf('%-18s %-14.6e %-14.6e %-12.6f %-12s %-12s\n','Euler explicito',E_f_exp,E_malla_exp,tiempo_exp,'no aplica','no aplica');

% Perfiles
tiempos_perfiles = [0,1,5,10];

figure;

for k = 1:length(tiempos_perfiles)

    tk = tiempos_perfiles(k);

    indice = round(tk/dt) + 1;

    subplot(2,2,k);

    plot(x_exp,T_exp(indice,:),'o-');

    hold on;

    plot(x_exp,T_exacta(x_exp,tk),'-','LineWidth',1.5);

    xlabel('x (m)');
    ylabel('T (C)');
    title(sprintf('t = %.1f s',tk));

    legend('Euler explicito','Solucion exacta','location','best');

    grid on;

end

print('perfiles_euler_explicito.png','-dpng');

% Estabilidad de Euler explicito
M_est = 50;
dx_est = L/M_est;
tf_est = 9.6;

mu_1 = 0.4;
mu_2 = 0.6;

dt_1 = mu_1*dx_est^2/alpha;
dt_2 = mu_2*dx_est^2/alpha;

% Condicion inicial perturbada
u0_pert = @(x) 80*sin(pi*x/L) + 1e-6*(-1).^round(x/dx_est);

% Caso mu = 0.4
[x_04,t_04,U_04] = euler_explicito(L,alpha,u0_pert,M_est,dt_1,tf_est);

T_04 = U_04 + 20;

% Caso mu = 0.6
[x_06,t_06,U_06] = euler_explicito(L,alpha,u0_pert,M_est,dt_2,tf_est);

T_06 = U_06 + 20;

% max |T_i^n - 20|
max_04 = max(abs(T_04 - 20),[],2);
max_06 = max(abs(T_06 - 20),[],2);

% Verificar caso inestable
indice_inestable = find(max_06 > 1e6,1);

if ~isempty(indice_inestable)

    tiempo_inestable = t_06(indice_inestable);

    fprintf('\nmu = 0.6 supera 10^6 C en t = %.6f s\n',tiempo_inestable);

end

% Grafica de estabilidad
figure;

semilogy(t_04,max_04,'-','LineWidth',1.5);

hold on;

semilogy(t_06,max_06,'-','LineWidth',1.5);

xlabel('Tiempo (s)');
ylabel('max |T_i^n - 20|');
title('Estabilidad de Euler explicito');

legend('mu = 0.4','mu = 0.6','location','best');

grid on;

print('estabilidad_euler_explicito.png','-dpng');

% Resumen de estabilidad
fprintf('\nEstabilidad de Euler explicito\n');
fprintf('mu = 0.4, dt = %.6f s: cumple mu <= 0.5\n',dt_1);
fprintf('mu = 0.6, dt = %.6f s: no cumple mu <= 0.5\n',dt_2);
