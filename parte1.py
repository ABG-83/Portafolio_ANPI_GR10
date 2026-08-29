import sympy as sp
import numpy as np


def secante(f, x0, x1, iterMax, tol):

    # Convertir la función ingresada como texto en una función numérica.
    x = sp.symbols("x")
    funcion_simbolica = sp.sympify(f)
    funcion_numerica = sp.lambdify(x, funcion_simbolica, "numpy")

    # Verificar si alguno de los valores iniciales es una raíz.
    if funcion_numerica(x0) == 0:
        xk = x0
        erk = 0
        k = 0
        conv = 1
        return xk, erk, k, conv

    if funcion_numerica(x1) == 0:
        xk = x1
        erk = 0
        k = 0
        conv = 1
        return xk, erk, k, conv

    # Inicializar las aproximaciones.
    x_anterior = x0
    xk = x1
    k = 0

    # Calcular el error inicial.
    erk = abs(funcion_numerica(xk))

    # Proceso iterativo.
    while k < iterMax and erk > tol:

        # Calcular el denominador del método de la secante.
        denominador = funcion_numerica(xk) - funcion_numerica(x_anterior)

        # Verificar que el denominador sea diferente de cero.
        if denominador == 0:
            conv = 0
            return xk, erk, k, conv

        # Incrementar el contador de iteraciones.
        k = k + 1

        # Calcular la nueva aproximación.
        x_nuevo = xk - funcion_numerica(xk)*(xk-x_anterior)/denominador

        # Actualizar las aproximaciones.
        x_anterior = xk
        xk = x_nuevo

        # Calcular el error.
        erk = abs(funcion_numerica(xk))

    # Determinar si el método convergió.
    if erk <= tol and k < iterMax:
        conv = 1
    else:
        conv = 0

    return xk, erk, k, conv

def biseccion(f, a, b, iterMax, tol):

    # Convertir la función ingresada como texto en una función numérica.
    x = sp.symbols("x")
    funcion_simbolica = sp.sympify(f)
    funcion_numerica = sp.lambdify(x, funcion_simbolica, "numpy")

    # Evaluar la función en los extremos del intervalo.
    valor_a = funcion_numerica(a)
    valor_b = funcion_numerica(b)

    # Verificar si alguno de los extremos es una raíz.
    if valor_a == 0:
        xk = a
        erk = 0
        k = 0
        conv = 1
        return xk, erk, k, conv

    if valor_b == 0:
        xk = b
        erk = 0
        k = 0
        conv = 1
        return xk, erk, k, conv

    # Verificar la condición necesaria del método.
    if valor_a*valor_b > 0:
        xk = np.nan
        erk = np.nan
        k = 0
        conv = 0
        return xk, erk, k, conv

    # Inicializar las variables.
    k = 0
    erk = np.inf
    xk = (a+b)/2

    # Proceso iterativo.
    while k < iterMax and erk > tol:

        # Incrementar el contador de iteraciones.
        k = k + 1

        # Calcular el punto medio.
        xk = (a+b)/2

        # Calcular el error.
        erk = abs(funcion_numerica(xk))

        # Verificar si se alcanzó la tolerancia.
        if erk <= tol:
            break

        # Seleccionar el nuevo intervalo.
        if funcion_numerica(a)*funcion_numerica(xk) < 0:
            b = xk
        else:
            a = xk

    # Determinar si el método convergió.
    if erk <= tol and k < iterMax:
        conv = 1
    else:
        conv = 0

    return xk, erk, k, conv

print(secante("x**2-2", 1, 2, 1000, 1e-8))
print(biseccion("x**2-2", 1, 2, 1000, 1e-8))

# Caso donde bisección no cumple cambio de signo
print(biseccion("x**2+1", -1, 1, 1000, 1e-8))

# Caso donde un valor inicial ya es raíz
print(secante("x**2-4", 2, 3, 1000, 1e-8))