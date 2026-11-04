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


% Persona 2: separar la factorizacion permite reutilizarla en todos los pasos.
function [multiplicadores, diagonal, superior] = factorizar_thomas(A)

    % Extraer las tres diagonales evita una eliminacion sobre la matriz completa.
    diagonal = diag(A);
    inferior = diag(A,-1);
    superior = diag(A,1);
    n = length(diagonal);
    multiplicadores = zeros(n-1,1);

    % La dominancia diagonal estricta garantiza pivotes no nulos en estas matrices.
    for i = 2:n
        multiplicadores(i-1) = inferior(i-1)/diagonal(i-1);
        diagonal(i) = diagonal(i) - multiplicadores(i-1)*superior(i-1);
    end

end


function z = metodo_thomas(multiplicadores, diagonal, superior, b)

    % Resolver el sistema triangular inferior usando los factores ya calculados.
    n = length(diagonal);
    y = b(:);
    for i = 2:n
        y(i) = y(i) - multiplicadores(i-1)*y(i-1);
    end

    % Recuperar la solucion por sustitucion hacia atras en el sistema superior.
    z = zeros(n,1);
    z(n) = y(n)/diagonal(n);
    for i = n-1:-1:1
        z(i) = (y(i) - superior(i)*z(i+1))/diagonal(i);
    end

end


function [x, t, U, residuos, conv] = euler_implicito(L, alpha, u0, M, dt, tf)

    % Preparar las mallas y el parametro que relaciona ambos pasos de discretizacion.
    dx = L/M;
    N = round(tf/dt);
    mu = alpha*dt/dx^2;
    x = (0:M)*dx;
    t = (0:N)*dt;
    tol = 1e-10;

    % Euler implicito resuelve A*u_nuevo = u_anterior, con A fija y tridiagonal.
    A = (1 + 2*mu)*eye(M-1);
    if M > 2
        secundaria = -mu*ones(M-2,1);
        A = A + diag(secundaria,1) + diag(secundaria,-1);
    end
    [multiplicadores, diagonal, superior] = factorizar_thomas(A);

    % Incluir los extremos nulos de u permite recuperar despues T = U + 20.
    u = u0(x(2:M));
    u = u(:);
    U = zeros(N+1,M+1);
    U(1,2:M) = u';
    residuos = zeros(N,1);
    conv = 1;

    % Reutilizar los factores y verificar cada sistema antes de aceptar el avance.
    for n = 1:N
        b = u;
        u_nuevo = metodo_thomas(multiplicadores,diagonal,superior,b);
        residuos(n) = norm(A*u_nuevo-b,2)/max(1,norm(b,2));
        if ~(residuos(n) < tol)
            conv = 0;
            warning('Euler implicito: fallo del sistema en el paso %d (t = %.6f s).',n,t(n+1));
            U = U(1:n,:);
            t = t(1:n);
            residuos = residuos(1:n);
            return;
        end
        u = u_nuevo;
        U(n+1,2:M) = u';
    end

end
