1;

function [x, t, U] = euler_explicito(L, alpha, u0, M, dt, tf)

    dx = L / M;
    N = round(tf / dt);
    mu = alpha * dt / dx^2;

    x = 0:dx:L;
    t = 0:dt:tf;

    U = zeros(N + 1, M + 1);

    % Condición inicial
    U(1, :) = u0(x);

    % Condiciones de frontera
    U(:, 1) = 0;
    U(:, M + 1) = 0;

    for n = 1:N
        for i = 2:M

            U(n + 1, i) = U(n, i) ...
                + mu * ( ...
                    U(n, i + 1) ...
                    - 2 * U(n, i) ...
                    + U(n, i - 1) ...
                );

        end
    end

end
