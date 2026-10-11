import time
import numpy as np


def metodo_qr(A, iterMax, tol):
    Ak = np.array(A, dtype=float, copy=True)
    n = Ak.shape[0]
    Vk = np.eye(n)
    I = np.eye(n)

    # Desplazamiento c
    C = np.linalg.norm(Ak, np.inf) + 1.0

    # Diagonal previa ordenada para medir el error
    dk = np.sort(np.diag(Ak))[::-1]

    con = 0
    k = 0

    for k in range(1, iterMax + 1):

        # Factorización QR de (Ak + C*I)
        Qk, Rk = np.linalg.qr(Ak + C * I)

        # Actualización de Ak y Vk
        Ak = Rk @ Qk - C * I
        Vk = Vk @ Qk

        # Nueva diagonal ordenada de la iteración actual
        dkN = np.sort(np.diag(Ak))[::-1]

        # Error como norma 2 entre las diagonales
        erk = np.linalg.norm(dkN - dk, ord=2)

        # Actualizar la diagonal previa
        dk = dkN

        # Criterio de parada
        if erk < tol:
            con = 1
            break

    # Reordenamiento de valores propios y columnas de vectores propios (de mayor a menor)
    lambdas = np.diag(Ak)
    idx = np.argsort(lambdas)[::-1]

    lambdas = lambdas[idx]
    Vk = Vk[:, idx]

    # Normalización de vectores propios 
    normas = np.linalg.norm(Vk, axis=0)
    normas[normas == 0] = 1.0
    Vk = Vk / normas

    return lambdas, Vk, erk, k, con


def metodo_jacobi(A, iterMax, tol, mostrar_progreso=False):

    # Copia de la matriz para no modificar la matriz original
    Ak = np.array(A, dtype=float, copy=True)

    # Tamaño de la matriz
    n = Ak.shape[0]

    # Matriz donde se acumularán los vectores propios
    Vk = np.eye(n)

    # Primera aproximación de los valores propios:
    # la diagonal de la matriz
    lambdas_prev = np.sort(np.diag(Ak))[::-1]

    # Inicialización
    erk = np.inf
    conv = 0
    k = 0

    # Medir el progreso sin cambiar las rotaciones ni la parada.
    total_pares = n * (n - 1) // 2
    inicio_metodo = time.perf_counter()
    for k in range(1, iterMax + 1):
        inicio_barrido = time.perf_counter()
        siguiente_aviso = 0.1
        if mostrar_progreso:
            error_previo = f"{erk:.6e}" if k > 1 else "aun no calculado"
            print(f"Barrido {k}: inicio; error previo {error_previo}; tol {tol:.1e}", flush=True)

        for p in range(n - 1):

            for q in range(p + 1, n):

                # Elemento fuera de la diagonal que se quiere eliminar
                apq = Ak[p, q]

                if abs(apq) <= 1e-15:
                    continue

                app = Ak[p, p]
                aqq = Ak[q, q]

                # Ángulo de la rotación de Jacobi
                theta = 0.5 * np.arctan2(
                    2.0 * apq,
                    app - aqq
                )

                c = np.cos(theta)
                s = np.sin(theta)

                # Rotación de las columnas p y q
                columna_p = Ak[:, p].copy()
                columna_q = Ak[:, q].copy()

                Ak[:, p] = c * columna_p + s * columna_q
                Ak[:, q] = -s * columna_p + c * columna_q

                # Rotación de las filas p y q
                fila_p = Ak[p, :].copy()
                fila_q = Ak[q, :].copy()

                Ak[p, :] = c * fila_p + s * fila_q
                Ak[q, :] = -s * fila_p + c * fila_q

                vector_p = Vk[:, p].copy()
                vector_q = Vk[:, q].copy()

                Vk[:, p] = c * vector_p + s * vector_q
                Vk[:, q] = -s * vector_p + c * vector_q


            # La estimacion corresponde al barrido actual, no al tiempo hasta converger.
            if mostrar_progreso:
                pares_recorridos = (p + 1) * (2 * n - p - 2) // 2
                fraccion = pares_recorridos / total_pares
                if fraccion >= siguiente_aviso or p == n - 2:
                    transcurrido = time.perf_counter() - inicio_barrido
                    restante = transcurrido * (1 - fraccion) / fraccion
                    print(f"Barrido {k}: {fraccion:.1%} de pares recorridos; "
                          f"{transcurrido:.1f} s; restante estimado del barrido "
                          f"{restante:.1f} s", flush=True)
                    siguiente_aviso += 0.1

        # Aproximaciones de los valores propios

        lambdas_nuevos = np.sort(np.diag(Ak))[::-1]

        erk = np.linalg.norm(
            lambdas_nuevos - lambdas_prev,
            ord=2
        )

        lambdas_prev = lambdas_nuevos
        # Mostrar el error permite compararlo con la tolerancia exigida.
        if mostrar_progreso:
            print(f"Barrido {k} terminado: error {erk:.6e}; tol {tol:.1e}; "
                  f"tiempo total {time.perf_counter()-inicio_metodo:.1f} s", flush=True)

        # Criterio de parada
        if erk < tol:
            conv = 1
            break


    # Ordenar valores propios y vectores propios conjuntamente

    diagonal_final = np.diag(Ak)

    indices = np.argsort(diagonal_final)[::-1]

    lambdas = diagonal_final[indices]
    Vk = Vk[:, indices]

    normas = np.linalg.norm(Vk, axis=0)

    normas[normas == 0] = 1.0

    Vk = Vk / normas

    return lambdas, Vk, erk, k, conv
    
