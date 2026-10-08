

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

% Comparar residuos reales con la escala de b y la tolerancia solicitada

fprintf("Los residuos absolutos estan entre %.8e y %.8e.\n",min(errores),max(errores));
residuos_relativos=errores/norm(b,2);
fprintf("Respecto a ||b||_2, estan entre %.8e y %.8e.\n",min(residuos_relativos),max(residuos_relativos));
disp("Estos valores permiten comparar que tan bien satisface cada aproximacion el sistema discreto.")
disp("Un residuo pequeno no significa una solucion algebraicamente exacta ni mide el error frente a la solucion continua.")
disp("La tolerancia es un criterio de parada de los iterativos; los directos no se detienen mediante ese criterio.")

% Contrastar los indicadores retornados con los residuos reales, sin cambiarlos

residuos_iter=[error_jacobi error_gauss_seidel error_cg];
convergencias_iter=[conv_jacobi conv_gauss_seidel conv_cg];
for i=1:3
  cumple_tol=residuos_iter(i)<tol;
  fprintf("%s: residuo=%.8e, k=%d, conv=%d, cumple residuo < tol: %d.\n",metodos_iter{i},residuos_iter(i),iteraciones_iter(i),convergencias_iter(i),cumple_tol);
  if convergencias_iter(i)~=cumple_tol
    disp("El indicador retornado no coincide con el residuo real; requiere revisar ese metodo.")
  elseif ~cumple_tol && iteraciones_iter(i)==iterMax
    disp("Se alcanzo el maximo de iteraciones sin satisfacer la tolerancia absoluta.")
  endif
endfor

disp(" ")
disp("- METODOS CON LOS MENORES TIEMPOS DE EJECUCION:")

% Obtener los menores tiempos de esta ejecucion dentro de cada grupo

[menor_directo,indice_directo]=min(tiempos(1:5));
[menor_iterativo,indice_iterativo]=min(tiempos(6:8));
fprintf(" - Menor tiempo entre los directos: %s con %.6f segundos.\n",metodos{indice_directo},menor_directo);
fprintf(" - Menor tiempo entre los iterativos: %s con %.6f segundos.\n",metodos_iter{indice_iterativo},menor_iterativo);
fprintf(" - Cholesky: %.6f s; LU: %.6f s; Gauss: %.6f s.\n",tiempo_cholesky,tiempo_lu,tiempo_gauss);
disp("El menor tiempo debe interpretarse junto con el residuo real y el cumplimiento de la tolerancia.")
disp("Los tiempos dependen del algoritmo, de los bucles o la vectorizacion y del equipo utilizado.")

disp(" ")
disp("- INFLUENCIA DE LA ESTRUCTURA TRIDIAGONAL SOBRE EL METODO DE THOMAS:")
disp("Thomas utiliza las tres diagonales de A y sus recurrencias requieren O(N) operaciones.")
disp("Gauss y LU, en las implementaciones generales utilizadas, recorren la matriz con costo O(N^3).")
disp("La estructura tridiagonal favorece a Thomas; la matriz A del programa se almacena como matriz llena.")

disp(" ")
disp("- VENTAJAS DE CHOLESKY BAJO LAS CONDICIONES DE LA MATRIZ:")
disp("A es simetrica, tiene diagonal positiva y es estrictamente dominante por filas: 510000 > 500000.")
disp("Estas propiedades garantizan que sea definida positiva y permiten aplicar Cholesky.")
disp("Cholesky obtiene A=L*L' calculando un solo factor triangular; LU calcula dos factores.")
disp("Su factorizacion requiere aproximadamente la mitad de operaciones que LU para matrices densas.")
disp("Esto no implica que el tiempo medido sea exactamente la mitad: aqui Cholesky esta vectorizado y LU usa bucles.")

disp(" ")
disp("- DIFERENCIAS OBSERVADAS ENTRE JACOBI Y GAUSS-SEIDEL:")

% Comparar iteraciones y tiempos sin deducir una velocidad universal

fprintf("Jacobi: %d iteraciones, %.6f s, residuo %.8e.\n",k_jacobi,tiempo_jacobi,error_jacobi);
fprintf("Gauss-Seidel: %d iteraciones, %.6f s, residuo %.8e.\n",k_gauss_seidel,tiempo_gauss_seidel,error_gauss_seidel);
disp("Jacobi utiliza la aproximacion anterior; Gauss-Seidel aprovecha las componentes nuevas mediante sustitucion hacia adelante.")
disp("El costo por iteracion es distinto: menos iteraciones no garantiza menor tiempo de ejecucion.")
disp("La dominancia diagonal garantiza convergencia en aritmetica exacta, pero el redondeo puede limitar el residuo alcanzable.")
disp("Llegar a iterMax no demuestra divergencia ni permite atribuir el resultado solamente a convergencia lenta.")

disp(" ")
disp("- AJUSTE DE LOS PARES ORDENADOS CON LA SOLUCION EXACTA:")
disp("La grafica permite corroborar visualmente los 501 puntos de diferencias finitas frente a la curva exacta.")
disp("En la grafica revisada, los puntos siguen muy de cerca la curva a la escala mostrada.")
disp("La superposicion visual no demuestra igualdad exacta: existen errores de discretizacion y de redondeo.")
disp("Las fronteras se fijan en 20 C y 50 C; las 499 temperaturas interiores se obtienen resolviendo AT=b.")
