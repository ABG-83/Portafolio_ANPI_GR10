1;
% Metodo de Eliminacion Gaussiana

function x=sol_elim_gauss(A,b)

  n=size(A,1);
  At=A;
  bt=b;

  % Eliminacion Gaussiana

  for k=1:n-1
    for i=k+1:n

      m=At(i,k)/At(k,k);

      for j=k:n
        At(i,j)=At(i,j)-m*At(k,j);
      endfor

      bt(i)=bt(i)-m*bt(k);

    endfor
  endfor

  % Sustitución hacia atrás

  x=zeros(n,1);

  for i=n:-1:1

    aux=0;

    for j=i+1:n
      aux=aux+At(i,j)*x(j);
    endfor

    x(i)=(1/At(i,i))*(bt(i)-aux);

  endfor

endfunction


% Método de Thomas

function x=metodo_thomas(A,d)

  b=diag(A);
  c=diag(A,1);
  a=diag(A,-1);
  a=([0 a'])';

  n=size(A,1);

  % Calcular vectores p y q
  p=zeros(n-1,1);
  q=zeros(n,1);

  p(1)=c(1)/b(1);
  q(1)=d(1)/b(1);

  for i=2:n-1

    aux=b(i)-p(i-1)*a(i);

    p(i)=c(i)/aux;
    q(i)=(d(i)-q(i-1)*a(i))/aux;

  endfor

  q(n)=(d(n)-q(n-1)*a(n))/(b(n)-p(n-1)*a(n));

  % Calcular vector x
  x=zeros(n,1);
  x(n)=q(n);

  for i=n-1:-1:1
    x(i)=q(i)-p(i)*x(i+1);
  endfor

endfunction


% Método de Jacobi

function [xk,erk,k,conv]=metodo_Jacobi(A,b,x0,tol,iterMax)

  xk=x0;

  d1=diag(A);
  D=diag(d1);
  Dinv=diag(1./d1);
  LmU=A-D;

  cj=Dinv*b;
  Tj=-Dinv*LmU;

  k=0;
  erk=norm(A*xk-b,2);
  conv=0;

  while erk>=tol && k<iterMax

    xk=Tj*xk+cj;

    erk=norm(A*xk-b,2);

    k=k+1;

  endwhile

  if erk<tol
    conv=1;
  endif

endfunction


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


function [x, Q, R] = metodo_QR(A, b)

  % tamaño del sistema
  n = size(A, 1);

  % matrices Q y R para rellenar
  Q = zeros(n, n);
  R = zeros(n, n);


  % Factorización QR
  for j = 1:n

    u = A(:, j); % Seleccionamos la columna j de la matriz A

    % Restamos la proyeccion de a_k sobre los q_j previos
    % El for no se ejecuta para j = 1, es decir u1 = u = A(:, j)
    for i = 1:j-1
      % Producto punto entre la columna original de A y la columna i de Q
      producto_punto = dot(A(:, j), Q(:, i));
      u = u - producto_punto * Q(:, i);
    end

    % Normalizamos el vector u para obtener q_j y agregarlo a Q
    Q(:, j) = u / norm(u);
  end


  % Calcular la matriz R
  R = Q' * A;


  % Resolver Ax = b
  % Calculamos c = Q^T * b
  c = Q' * b;

  % Resolvemos el sistema Rx = c utilizando sustitución hacia atrás
  x=sust_atras(R,c);

end



function [xk, erk, k, conv] = gradiente_conjugado(A, b, x0, iterMax, tol)

  xk = x0;
  rk = b - A * xk;   % residual inicial: cuánto le falta a Ax para llegar a b
  pk = rk;           % la primera dirección es el propio residual

  k = 0;
  erk = norm(rk, 2); % error inicial medido con norma 2
  conv = 0;

  while k < iterMax

    Apk = A * pk;  % se guarda A*pk porque se usa dos veces

    % el tamaño de cuánto se debe avanzar en la dirección pk
    alpha = (rk' * rk) / (pk' * Apk);

    xk = xk + alpha * pk;  % nueva aproximación

    rk_nuevo = rk - alpha * Apk;

    k = k + 1;

    erk = norm(rk_nuevo, 2);

    if erk < tol
      conv = 1;
      rk = rk_nuevo;
      break;
    end

    % beta ajusta la nueva dirección para no repetir direcciones anteriores
    beta = (rk_nuevo' * rk_nuevo) / (rk' * rk);

    % nueva dirección: residual + un poco de la dirección anterior
    pk = rk_nuevo + beta * pk;

    rk = rk_nuevo;  % se actualiza el residual para la siguiente vuelta

  end

end

