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
