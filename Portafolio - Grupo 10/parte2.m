

clc;
clear;
close all;

source("parte1.m");

% Parte II - Comparacion computacional de los metodos


% Tamano del sistema

n=300;


% Matriz A

diagonal_principal=4*ones(n,1);
diagonal_secundaria=-1*ones(n-1,1);

A=diag(diagonal_principal)+diag(diagonal_secundaria,1)+diag(diagonal_secundaria,-1);


% Vector b

b=ones(n,1);


% Parametros de los metodos iterativos

x0=zeros(n,1);
iterMax=10000;
tol=1e-8;


% Eliminacion Gaussiana

tic;
x_gauss=sol_elim_gauss(A,b);
tiempo_gauss=toc;

error_gauss=norm(A*x_gauss-b,2);


% Metodo de Thomas

tic;
x_thomas=metodo_thomas(A,b);
tiempo_thomas=toc;

error_thomas=norm(A*x_thomas-b,2);


% Metodo de Jacobi

tic;
[x_jacobi,erk_jacobi,k_jacobi,conv_jacobi]=metodo_Jacobi(A,b,x0,tol,iterMax);
tiempo_jacobi=toc;

error_jacobi=norm(A*x_jacobi-b,2);


% Factorizacion LU

tic;
x_lu=sol_LU(A,b);
tiempo_lu=toc;

error_lu=norm(A*x_lu-b,2);


% Cholesky

tic;
x_cholesky=sol_Cholesky(A,b);
tiempo_cholesky=toc;

error_cholesky=norm(A*x_cholesky-b,2);


% QR
tic;
[x_qr, Q, R] = metodo_QR(A, b);
tiempo_qr= toc;
error_qr = norm(A * x_qr - b, 2);


% Gauss-Seidel

tic;
[x_gauss_seidel,erk_gauss_seidel,k_gauss_seidel,conv_gauss_seidel]=metodo_Gauss_Seidel(A,b,x0,tol,iterMax);
tiempo_gauss_seidel=toc;

error_gauss_seidel=norm(A*x_gauss_seidel-b,2);


% Gradiente Conjugado

tic;
[x_cg, erk_cg, k_cg, conv_cg] = gradiente_conjugado(A, b, x0, iterMax, tol);
tiempo_grad_conj = toc;
error_grad_conj= norm(A * x_cg - b, 2);





% -------- Tabla Comparativa Parte II --------

disp(" ")
disp("Tabla Comparativa de Metodos")
fprintf("%-22s %18s %16s %12s %8s\n","Metodo","Residuo ||Ax-b||_2","Tiempo (s)","Iteraciones","conv");
fprintf("%-22s %18.8e %16.6f %12s %8s\n","Eliminacion Gaussiana",error_gauss,tiempo_gauss,"-","-");
fprintf("%-22s %18.8e %16.6f %12s %8s\n","Factorizacion LU",error_lu,tiempo_lu,"-","-");
fprintf("%-22s %18.8e %16.6f %12s %8s\n","Cholesky",error_cholesky,tiempo_cholesky,"-","-");
fprintf("%-22s %18.8e %16.6f %12s %8s\n","QR",error_qr,tiempo_qr,"-","-");
fprintf("%-22s %18.8e %16.6f %12s %8s\n","Thomas",error_thomas,tiempo_thomas,"-","-");
fprintf("%-22s %18.8e %16.6f %12d %8d\n","Jacobi",error_jacobi,tiempo_jacobi,k_jacobi,conv_jacobi);
fprintf("%-22s %18.8e %16.6f %12d %8d\n","Gauss-Seidel",error_gauss_seidel,tiempo_gauss_seidel,k_gauss_seidel,conv_gauss_seidel);
fprintf("%-22s %18.8e %16.6f %12d %8d\n","Gradiente Conjugado",error_grad_conj,tiempo_grad_conj,k_cg,conv_cg);


% -------- Graficas Comparativas --------

metodos={"Elim. Gaussiana","LU","Cholesky","QR","Thomas","Jacobi","Gauss-Seidel","Grad. Conjugado"};
errores=[error_gauss error_lu error_cholesky error_qr error_thomas error_jacobi error_gauss_seidel error_grad_conj];
tiempos=[tiempo_gauss tiempo_lu tiempo_cholesky tiempo_qr tiempo_thomas tiempo_jacobi tiempo_gauss_seidel tiempo_grad_conj];

% 1. Comparacion de Errores (Escala logaritmica)
figure;
semilogy(1:8,errores,"o-","LineWidth",2,"MarkerSize",6);
set(gca,"XTick",1:8,"XTickLabel",metodos);
grid on;
xlabel("Metodo");
ylabel("Residuo ||Ax-b||_2");
title("Comparacion de Errores");

% 2. Comparacion de Tiempos de Ejecucion (Escala logaritmica)
figure;
semilogy(1:8,tiempos,"s-","LineWidth",2,"MarkerSize",6);
set(gca,"XTick",1:8,"XTickLabel",metodos);
grid on;
xlabel("Metodo");
ylabel("Tiempo de ejecucion (s)");
title("Comparacion de Tiempos de Ejecucion");

% 3. Comparacion de Iteraciones (Metodos Iterativos)
metodos_iter={"Jacobi","Gauss-Seidel","Grad. Conjugado"};
iteraciones_iter=[k_jacobi k_gauss_seidel k_cg];

figure;
bar(iteraciones_iter,0.5);
set(gca,"XTick",1:3,"XTickLabel",metodos_iter);
grid on;
xlabel("Metodo Iterativo");
ylabel("Numero de Iteraciones (k)");
title("Iteraciones Requeridas por Metodos Iterativos");


% -------- Analisis Comparativo --------

disp(" ")
disp("================== ANALISIS COMPARATIVO - PARTE II =======================")
disp("Estructura de la Matriz A:")
disp("  La matriz A es tridiagonal, simetrica y estrictamente dominante por diagonal")
disp("  (|4| > |-1| + |-1|). Esto garantiza la convergencia de los metodos de Jacobi")
disp("  y Gauss-Seidel, y asegura que A sea definida positiva, permitiendo el uso de")
disp("  Cholesky y Gradiente Conjugado.")
disp(" ")
disp("Tiempos de Ejecucion:")
disp("- Thomas es por amplio margen el metodo mas rapido, ya que explota la estructura")
disp("  tridiagonal reduciendo la complejidad computacional a O(N).")
disp("- Cholesky es aproximadamente el doble de rapido que LU y Eliminacion Gaussiana,")
disp("  ya que solo calcula una matriz L y aprovecha la simetria de A.")
disp("- QR es el mas lento de los directos debido al costo de Gram-Schmidt sobre matrices llenas.")
disp(" ")
disp("Metodos Iterativos:")
disp("- Gauss-Seidel converge sustancialmente mas rapido que Jacobi (requiere cerca de ")
disp("  la mitad de iteraciones) porque utiliza inmediatamente los valores actualizados.")
disp("- Gradiente Conjugado destaca por requerir un numero reducido de iteraciones gracias")
disp("  a la busqueda en direcciones A-conjugadas para matrices simetricas definidas positivas.")

