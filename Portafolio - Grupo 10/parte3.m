

clc;
clear;
close all;

source("parte1.m");


% -------- Construccion del sistema de temperaturas --------

% Puntos de la barra y condiciones de frontera

h=1/500;
x=(0:500)'*h;
x_interior=x(2:500);
n=499;

T0=20;
T500=50;


% Matriz A de las 499 temperaturas interiores

diagonal_principal=510000*ones(n,1);
diagonal_secundaria=-250000*ones(n-1,1);

A=diag(diagonal_principal)+diag(diagonal_secundaria,1)+diag(diagonal_secundaria,-1);


% Vector b con la fuente de calor y los aportes de las fronteras

g=300000*x_interior+(100000+10*pi^2)*sin(pi*x_interior)+200000;
b=g;
b(1)=b(1)+250000*T0;
b(n)=b(n)+250000*T500;


% Parametros del metodo iterativo

T_inicial=zeros(n,1);
iterMax=10000;
tol=1e-8;


% -------- Metodos de Persona 1 --------


% Eliminacion Gaussiana

tic;
T_gauss=sol_elim_gauss(A,b);
tiempo_gauss=toc;

error_gauss=norm(A*T_gauss-b,2);

temperaturas_gauss=[T0;T_gauss;T500];
pares_gauss=[x temperaturas_gauss];


% Metodo de Thomas

tic;
T_thomas=metodo_thomas(A,b);
tiempo_thomas=toc;

error_thomas=norm(A*T_thomas-b,2);

temperaturas_thomas=[T0;T_thomas;T500];
pares_thomas=[x temperaturas_thomas];


% Metodo de Jacobi

tic;
[T_jacobi,erk_jacobi,k_jacobi,conv_jacobi]=metodo_Jacobi(A,b,T_inicial,tol,iterMax);
tiempo_jacobi=toc;

error_jacobi=norm(A*T_jacobi-b,2);

temperaturas_jacobi=[T0;T_jacobi;T500];
pares_jacobi=[x temperaturas_jacobi];

% -------- Metodos de Persona 2 --------

% Factorizacion LU

tic;
T_lu=sol_LU(A,b);
tiempo_lu=toc;

error_lu=norm(A*T_lu-b,2);
temperaturas_lu=[T0;T_lu;T500];
pares_lu=[x temperaturas_lu];


% Cholesky

tic;
T_cholesky=sol_Cholesky(A,b);
tiempo_cholesky=toc;

error_cholesky=norm(A*T_cholesky-b,2);
temperaturas_cholesky=[T0;T_cholesky;T500];
pares_cholesky=[x temperaturas_cholesky];


% Gauss-Seidel

tic;
[T_gauss_seidel,erk_gauss_seidel,k_gauss_seidel,conv_gauss_seidel]=metodo_Gauss_Seidel(A,b,T_inicial,tol,iterMax);
tiempo_gauss_seidel=toc;

error_gauss_seidel=norm(A*T_gauss_seidel-b,2);
temperaturas_gauss_seidel=[T0;T_gauss_seidel;T500];
pares_gauss_seidel=[x temperaturas_gauss_seidel];


% -------- Resultados de los metodos integrados --------

% Tabla de errores, tiempos e iteraciones  Persona 2

disp(" ")
disp("Parte III - Resultados de Persona 2")
fprintf("%-16s %16s %16s %12s %8s\n","Metodo","Residuo norma 2","Tiempo (s)","Iteraciones","conv");
fprintf("%-16s %16.8e %16.6f %12s %8s\n","LU",error_lu,tiempo_lu,"-","-");
fprintf("%-16s %16.8e %16.6f %12s %8s\n","Cholesky",error_cholesky,tiempo_cholesky,"-","-");
fprintf("%-16s %16.8e %16.6f %12d %8d\n","Gauss-Seidel",error_gauss_seidel,tiempo_gauss_seidel,k_gauss_seidel,conv_gauss_seidel);


% Comparacion de residuos y tiempos de los tres metodos

metodos={"LU","Cholesky","Gauss-Seidel"};
errores=[error_lu error_cholesky error_gauss_seidel];
tiempos=[tiempo_lu tiempo_cholesky tiempo_gauss_seidel];

figure;
semilogy(1:3,errores,"o-","LineWidth",2);
set(gca,"XTick",1:3,"XTickLabel",metodos);
xlim([0.5 3.5]);
xlabel("Metodo");
ylabel("Residuo ||AT-b||_2");
title("Barra metalica - Residuos de Persona 2");
grid on;

figure;
semilogy(1:3,tiempos,"o-","LineWidth",2);
set(gca,"XTick",1:3,"XTickLabel",metodos);
xlim([0.5 3.5]);
xlabel("Metodo");
ylabel("Tiempo de ejecucion (s)");
title("Barra metalica - Tiempos de Persona 2");
grid on;


% Iteraciones del metodo iterativo integrado por Persona 2

figure;
bar(1,k_gauss_seidel);
set(gca,"XTick",1,"XTickLabel",{"Gauss-Seidel"});
xlim([0.5 1.5]);
xlabel("Metodo");
ylabel("Iteraciones ejecutadas");
title("Barra metalica - Iteraciones de Gauss-Seidel");
grid on;


% Corroboracion grafica con la solucion exacta y los 501 puntos de Cholesky

x_exacta=linspace(0,1,2001);
T_exacta=20+30*x_exacta+10*sin(pi*x_exacta);

figure;
plot(x_exacta,T_exacta,"b-","LineWidth",2);
hold on;
scatter(x,temperaturas_cholesky,10,"r");
xlabel("Posicion x (m)");
ylabel("Temperatura (C)");
title("Distribucion de temperatura en la barra metalica");
legend("Solucion exacta","Diferencias finitas - Cholesky","Location","southeast");
grid on;
hold off;


% -------- Analisis de los resultados de Persona 2 --------

% Relacionar las mediciones con las propiedades de los metodos y de A

disp(" ")
disp("Analisis de los metodos de Persona 2")
fprintf("Los residuos de LU, Cholesky y Gauss-Seidel son %.8e, %.8e y %.8e.\n",errores);
disp("Estos valores miden que tan bien se satisface AT=b; no son el error respecto a la temperatura exacta.")
[tiempo_menor,indice_menor]=min(tiempos);
fprintf("En esta ejecucion, el menor tiempo entre los tres metodos fue el de %s: %.6f s.\n",metodos{indice_menor},tiempo_menor);
disp("A es simetrica, tiene diagonal positiva y es estrictamente diagonalmente dominante: 510000 > 500000.")
disp("Por estas propiedades A es definida positiva y se puede aplicar Cholesky.")
disp("Cholesky calcula un solo factor triangular L y utiliza su transpuesta; LU construye L y U.")
disp("La version vectorizada de Cholesky tambien influye en los tiempos medidos frente a los bucles de LU.")
disp("Gauss-Seidel utiliza las componentes nuevas al resolver el sistema triangular en cada iteracion.")
fprintf("Gauss-Seidel ejecuto %d iteraciones; conv=%d.\n",k_gauss_seidel,conv_gauss_seidel);
if conv_gauss_seidel==1
  disp("Gauss-Seidel alcanzo el residuo menor que 1e-8 solicitado.")
else
  disp("Gauss-Seidel alcanzo el maximo de iteraciones sin obtener un residuo menor que 1e-8.")
  disp("La convergencia teorica no elimina el limite de precision de la aritmetica de punto flotante.")
endif
disp("La grafica superpone la curva exacta y los 501 puntos numericos para corroborar su comportamiento.")
disp("Los extremos se fijaron en 20 C y 50 C; las 499 temperaturas interiores provienen del sistema AT=b.")
disp("La comparacion completa con los otros cinco metodos queda pendiente de su integracion por el grupo.")
