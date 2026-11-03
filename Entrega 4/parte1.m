1;

function [x, t, U] = euler_explicito(L, alpha, u0, M, dt, tf)

    dx = L/M;
    N = round(tf/dt);
    mu = alpha*dt/dx^2;

    x = (0:M)*dx;
    t = (0:N)*dt;

    % Matriz A
    A = (1 - 2*mu)*eye(M - 1);

    for i = 1:M-2
        A(i,i+1) = mu;
        A(i+1,i) = mu;
    end

    % Condicion inicial
    u = u0(x(2:M));
    u = u(:);

    % Guardar solucion
    U = zeros(N + 1,M + 1);
    U(1,2:M) = u';

    if mu > 0.5
        warning('El metodo puede ser inestable');
    end

    % Iteracion temporal
    for n = 1:N
        u = A*u;
        U(n+1,2:M) = u';
    end

end
