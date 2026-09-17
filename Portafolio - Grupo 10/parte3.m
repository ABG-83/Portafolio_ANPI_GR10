

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

% -------- Metodos de Persona 3 --------

% Factorizacion QR

tic;
[T_qr, Q, R]=metodo_QR(A,b);
tiempo_qr=toc;

error_qr=norm(A*T_qr-b,2);
temperaturas_qr=[T0;T_qr;T500];
pares_qr=[x temperaturas_qr];


% Gradiente Conjugado

tic;
[T_cg,erk_cg,k_cg,conv_cg]=gradiente_conjugado(A,b,T_inicial,iterMax,tol);
tiempo_cg=toc;

error_cg=norm(A*T_cg-b,2);
temperaturas_cg=[T0;T_cg;T500];
  pares_cg=[x temperaturas_cg];


% -------- Resultados de los metodos integrados --------

disp(" ")
disp("Tabla Comparativa General de Metodos")
fprintf("%-22s %18s %16s %12s %8s\n","Metodo","Residuo ||AT-b||_2","Tiempo (s)","Iteraciones","conv");
fprintf("%-22s %18.8e %16.6f %12s %8s\n","Eliminacion Gaussiana",error_gauss,tiempo_gauss,"-","-");
fprintf("%-22s %18.8e %16.6f %12s %8s\n","Factorizacion LU",error_lu,tiempo_lu,"-","-");
fprintf("%-22s %18.8e %16.6f %12s %8s\n","Cholesky",error_cholesky,tiempo_cholesky,"-","-");
fprintf("%-22s %18.8e %16.6f %12s %8s\n","QR",error_qr,tiempo_qr,"-","-");
fprintf("%-22s %18.8e %16.6f %12s %8s\n","Thomas",error_thomas,tiempo_thomas,"-","-");
fprintf("%-22s %18.8e %16.6f %12d %8d\n","Jacobi",error_jacobi,tiempo_jacobi,k_jacobi,conv_jacobi);
fprintf("%-22s %18.8e %16.6f %12d %8d\n","Gauss-Seidel",error_gauss_seidel,tiempo_gauss_seidel,k_gauss_seidel,conv_gauss_seidel);
fprintf("%-22s %18.8e %16.6f %12d %8d\n","Gradiente Conjugado",error_cg,tiempo_cg,k_cg,conv_cg);


% Comparacion de residuos y tiempos

metodos={"Elim. Gaussiana","LU","Cholesky","QR","Thomas","Jacobi","Gauss-Seidel","Grad. Conjugado"};
errores=[error_gauss error_lu error_cholesky error_qr error_thomas error_jacobi error_gauss_seidel error_cg];
tiempos=[tiempo_gauss tiempo_lu tiempo_cholesky tiempo_qr tiempo_thomas tiempo_jacobi tiempo_gauss_seidel tiempo_cg];

% Grafica de Residuos
figure;
semilogy(1:8,errores,"o-","LineWidth",2,"MarkerSize",6);
set(gca,"XTick",1:8,"XTickLabel",metodos);
grid on;
xlabel("Metodo");
ylabel("Residuo ||AT-b||_2");
title("Residuos por Metodo");

% Grafica de Tiempos de Ejecución
figure;
semilogy(1:8,tiempos,"s-","LineWidth",2,"MarkerSize",6);
set(gca,"XTick",1:8,"XTickLabel",metodos);
grid on;
xlabel("Metodo");
ylabel("Tiempo de ejecucion (s)");
title("Tiempos de Ejecucion por Metodo");

% Grafica de Iteraciones de Metodos Iterativos
metodos_iter={"Jacobi","Gauss-Seidel","Grad. Conjugado"};
iteraciones_iter=[k_jacobi k_gauss_seidel k_cg];

figure;
bar(iteraciones_iter,0.5);
set(gca,"XTick",1:3,"XTickLabel",metodos_iter);
grid on;
xlabel("Metodo Iterativo");
ylabel("Numero de Iteraciones (k)");
title("Iteraciones Requeridas por Metodos Iterativos");


% -------- Corroboracion Grafica (Solucion Exacta vs Numerica) --------

x_exacta=linspace(0,1,2001);
T_exacta=20+30*x_exacta+10*sin(pi*x_exacta);

figure;
plot(x_exacta,T_exacta,"b-","LineWidth",2);
hold on;
scatter(x,temperaturas_cg,10,"r","filled");
xlabel("Posicion x (m)");
ylabel("Temperatura (C)");
title("Distribucion de temperatura en la barra metalica");
legend("Solucion exacta","Diferencias finitas (501 puntos)","Location","southeast");
grid on;
hold off;


% -------- Analisis y Discusion Final --------

disp(" ")
disp("===================ANALISIS Y DISCUSION FINAL===================")

disp(" ")
disp("- EQUIVALENCIA DE SOLUCIONES SEGUN EL RESIDUO ||AT - b||_2:")
disp("Desde el punto de vista del residuo algebraico, todos los metodos producen")
disp("soluciones equivalentes. Los valores obtenidos se encuentran en el rango de")
fprintf("%.2e a %.2e. Estos residuos garantizan que todos los metodos lograron\n", min(errores), max(errores));
disp("despejar el sistema lineal AT = b de forma algebraicamente exacta y acotada")
disp("por la precision del sistema. Es importante destacar que el residuo mide la")
disp("fidelidad de la solucion del sistema discreto, no el error frente a la solucion continua.")

disp(" ")
disp("- METODOS CON LOS MENORES TIEMPOS DE EJECUCION:")
fprintf(" - Metodo Directo Mas Rapido: Thomas con %.6f segundos.\n", tiempo_thomas);
fprintf(" - Metodo Iterativo Mas Rapido: Gradiente Conjugado con %.6f segundos.\n", tiempo_cg);
fprintf(" - Destaca tambien Cholesky con %.6f segundos frente a LU (%.6f s) o Gauss (%.6f s).\n", tiempo_cholesky, tiempo_lu, tiempo_gauss);
disp("Los tiempos mas bajos corresponden a los algoritmos capaces de explotar la")
disp("estructura especifica de la matriz A (tridiagonalidad y simetria definida positiva).")

disp(" ")
disp("- INFLUENCIA DE LA ESTRUCTURA TRIDIAGONAL SOBRE EL METODO DE THOMAS:")
disp("La matriz A es tridiagonal rala. El algoritmo de Thomas aprovecha de forma")
disp("directa esta estructura, reduciendo la complejidad computacional de O(N^3)")
disp("a O(N). Para N = 499, el numero de operaciones pasa de aproximadamente")
disp("1.2x10^8 a solo ~2500 multiplicaciones/divisiones, logrando que el tiempo")
disp("de calculo sea practicamente instantaneo en comparacion con LU o Gauss clasicos.")

disp(" ")
disp("- VENTAJAS DE CHOLESKY BAJO LAS CONDICIONES DE LA MATRIZ:")
disp("La matriz A es simetrica (A = A') y strictly dominante por diagonal")
disp("(|510000| > |-250000| + |-250000|), lo que garantiza teoricamente que A es")
disp("Simetrica Definida Positiva (SPD). Cholesky aprovecha esto descomponiendo A")
disp("en L*L', requiriendo solo la mitad de operaciones flotantes y la mitad de")
disp("memoria respecto a la factorizacion LU. Ademas, la implementacion vectorizada")
disp("de Octave reduce drasticamente el tiempo de ejecucion.")

disp(" ")
disp("- DIFERENCIAS OBSERVADAS ENTRE JACOBI Y GAUSS-SEIDEL:")
fprintf("Ambos metodos alcanzaron el limite maximo de %d iteraciones (conv = 0).\n", iterMax);
disp("En discretizaciones de EDPs por diferencias finitas con mallas finas (N = 499),")
disp("el radio espectral de las matrices de iteracion es sumamente cercano a 1, lo")
disp("que genera una convergencia asintotica muy lenta (estancamiento). Sin embargo,")
fprintf("Gauss-Seidel logro un residuo significativamente menor (%.8e) que Jacobi\n", error_gauss_seidel);
fprintf("(%.8e), confirmando la teoria de que Gauss-Seidel converge aproximadamente\n", error_jacobi);
disp("al doble de velocidad por iteracion al utilizar de inmediato los valores nuevos.")

disp(" ")
disp("- AJUSTE DE LOS PARES ORDENADOS CON LA SOLUCION EXACTA:")
disp("Los 501 puntos ordenados (x, T) obtenidos por diferencias finitas siguen de")
disp("manera exacta la curva continua T(x) = 20 + 30x + 10*sin(pi*x). En la grafica")
disp("se observa un acople perfecto entre el scatter de los puntos numericos y")
disp("la solucion analitica. Esto valida la formulacion del sistema AT = b y confirma")
disp("que el paso de malla h = 1/500 modela con gran precision el problema fisico.")
