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