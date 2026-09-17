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


% Cholesky


% QR


% Gauss-Seidel


% Gradiente Conjugado

