

clc;
clear;

source("parte1.m");


% Prueba Eliminacion Gaussiana

disp(" ")
disp("Prueba Eliminacion Gaussiana")

A=[10 -1 2;
   -1 11 -1;
    2 -1 10];

b=[6 22 -10].';

x_gauss=sol_elim_gauss(A,b)

error_gauss=norm(A*x_gauss-b)


% Prueba Thomas

disp(" ")
disp("Prueba Thomas")

A=[ 2 -1  0  0;
   -1  3 -1  0;
    0 -1  5 -1;
    0  0 -1  3];

d=[1 1 3 2].';

x_thomas=metodo_thomas(A,d)

error_thomas=norm(A*x_thomas-d)


% Prueba Jacobi

disp(" ")
disp("Prueba Jacobi")

A=[10 -1 2;
   -1 11 -1;
    2 -1 10];

b=[6 22 -10].';

x0=[0 0 0].';

tol=1e-10;
iterMax=1000;

[x_jacobi,erk,k,conv]=metodo_Jacobi(A,b,x0,tol,iterMax)




% Funciones correspondientes a Persona 2.

% -------- Factorizacion LU --------

function x=sol_LU(A,b)

  % Obtener L y U sin pivoteo, guardando los multiplicadores de Gauss

  n=size(A,1);
  U=A;
  L=eye(n);

  for k=1:n-1
    for i=k+1:n
      m=U(i,k)/U(k,k);
      L(i,k)=m;
      for j=k:n
        U(i,j)=U(i,j)-m*U(k,j);
      endfor
    endfor
  endfor

  % Resolver Ly=b y luego Ux=y

  y=sust_adelante(L,b);
  x=sust_atras(U,y);

endfunction


% -------- Factorizacion de Cholesky --------

function x=sol_Cholesky(A,b)

  % Construir L por columnas como en fact_cholesky_f del ejemplo de clase

  n=size(A,1);
  L=zeros(n,n);

  for i=1:n

    % Calcular el elemento diagonal con la suma de cuadrados anteriores

    L(i,i)=sqrt(A(i,i)-L(i,1:i-1)*L(i,1:i-1)');

    % Calcular los elementos debajo de la diagonal de la columna i

    if i<n
      L(i+1:n,i)=(A(i+1:n,i)-L(i+1:n,1:i-1)*L(i,1:i-1)')/L(i,i);
    endif
  endfor

  % Resolver Ly=b y luego L'x=y

  y=sust_adelante(L,b);
  x=sust_atras(L',y);

endfunction


% -------- Metodo de Gauss-Seidel --------

function [xk,erk,k,conv]=metodo_Gauss_Seidel(A,b,x0,tol,iterMax)

  % Separar A=L+D+U y formar la matriz triangular inferior M=L+D

  L=tril(A,-1);
  D=diag(diag(A));
  U=triu(A,1);
  M=L+D;

  % Inicializar la aproximacion y calcular la norma del residuo

  xk=x0;
  erk=norm(A*xk-b,2);
  k=0;

  % Resolver M*x_nuevo=b-U*xk por sustitucion hacia adelante

  while erk>=tol && k<iterMax
    c=b-U*xk;
    xk=sust_adelante(M,c);
    erk=norm(A*xk-b,2);
    k=k+1;
  endwhile

  % Indicar si se alcanzo la tolerancia estricta del enunciado

  conv=double(erk<tol);

endfunction


% -------- Sustituciones triangulares de Persona 2 --------

function y=sust_adelante(A,b)

  % Reservar el vector solucion del sistema triangular inferior

  n=size(A,1);
  y=zeros(n,1);

  % Calcular cada componente usando las componentes anteriores

  for i=1:n
    aux=0;
    if i>1
      aux=A(i,1:i-1)*y(1:i-1);
    endif
    y(i)=(b(i)-aux)/A(i,i);
  endfor

endfunction


function x=sust_atras(A,b)

  % Reservar el vector solucion del sistema triangular superior

  n=size(A,1);
  x=zeros(n,1);

  % Calcular cada componente desde la ultima fila hacia la primera

  for i=n:-1:1
    aux=0;
    if i<n
      aux=A(i,i+1:n)*x(i+1:n);
    endif
    x(i)=(b(i)-aux)/A(i,i);
  endfor

endfunction
