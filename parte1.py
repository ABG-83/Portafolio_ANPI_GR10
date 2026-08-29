"""Implementa metodos numericos para aproximar raices de funciones."""

import sympy as sp
import numpy as np


def secante(f, x0, x1, iterMax, tol):
    x = sp.symbols("x")
    funcion_simbolica = sp.sympify(f)
    funcion_numerica = sp.lambdify(x, funcion_simbolica, "numpy")

# se verifica si algun valor inciial es una raiz 
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

#aprox
    x_anterior = x0
    xk = x1
    k = 0

#error inicial
    erk = abs(funcion_numerica(xk))


    while k < iterMax and erk > tol:
        denominador = funcion_numerica(xk) - funcion_numerica(x_anterior)
        if denominador == 0:
            conv = 0
            return xk, erk, k, conv

        k = k + 1

#nueva aprox
        x_nuevo = xk - funcion_numerica(xk)*(xk-x_anterior)/denominador

#actualiza las aprox
        x_anterior = xk
        xk = x_nuevo

#calcular el error
        erk = abs(funcion_numerica(xk))

    if erk <= tol and k < iterMax:
        conv = 1
    else:
        conv = 0

    return xk, erk, k, conv



def biseccion(f, a, b, iterMax, tol):
    x = sp.symbols("x")
    funcion_simbolica = sp.sympify(f)
    funcion_numerica = sp.lambdify(x, funcion_simbolica, "numpy")

#evaluar en extremos
    valor_a = funcion_numerica(a)
    valor_b = funcion_numerica(b)

#verificar si algun extremo es raiz
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

    if valor_a*valor_b > 0:
        xk = np.nan
        erk = np.nan
        k = 0
        conv = 0
        return xk, erk, k, conv

    k = 0
    erk = np.inf
    xk = (a+b)/2

    while k < iterMax and erk > tol:
        k = k + 1
        xk = (a+b)/2
        erk = abs(funcion_numerica(xk))

#verifica si se alcanza la tol
        if erk <= tol:
            break

        if funcion_numerica(a)*funcion_numerica(xk) < 0:
            b = xk
        else:
            a = xk

#si el metodo converge 
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



# Metodo de Newton Raphson


def newton_raphson(f, x0, iterMax, tol):
    x = sp.symbols("x")
    funcion_simbolica = sp.sympify(f)
    derivada_simbolica = sp.diff(funcion_simbolica, x)
    funcion_numerica = sp.lambdify(x, funcion_simbolica, "numpy")
    derivada_numerica = sp.lambdify(x, derivada_simbolica, "numpy")

#aprox y error inicial
    xk = x0
    k = 0
    erk = abs(funcion_numerica(xk))

#se verifica si el valor inicial es una raiz
    if erk <= tol:
        conv = 1
        return xk, erk, k, conv

#calcular las aproximaciones
    while k < iterMax and erk > tol:
        valor_derivada = derivada_numerica(xk)

#se evita dividir entre cero
        if valor_derivada == 0:
            conv = 0
            return xk, erk, k, conv

        k = k + 1

#nueva aprox
        xk = xk - funcion_numerica(xk)/valor_derivada

#calcular el error
        erk = abs(funcion_numerica(xk))

#se verifica si el metodo converge
    if erk <= tol and k < iterMax:
        conv = 1
    else:
        conv = 0

    return xk, erk, k, conv


# Prueba 
print(newton_raphson("x**2-2", 1, 1000, 1e-8))

 
