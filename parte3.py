from parte1 import secante, biseccion, newton_raphson, steffensen, falsa_posicion, muller

import sympy as sp
import numpy as np
import time
import matplotlib.pyplot as plt


def analisis_funcion(g_str):
    x = sp.symbols("x")
    funcion_simbolica = sp.sympify(g_str)
    funcion_numerica = sp.lambdify(x, funcion_simbolica, "numpy")

    valores_f = np.linspace(0.005, 0.05, 1000)
    valores_g = funcion_numerica(valores_f)

    # Graficar la función
    plt.figure(figsize=(10, 6))
    plt.plot(valores_f, valores_g, label="g(f)")
    plt.axhline(0)
    plt.xlabel("Factor de fricción f")
    plt.ylabel("g(f)")
    plt.title("Análisis de la ecuación de Colebrook-White")
    plt.grid(True)
    plt.legend()
    plt.show()

    # Valores iniciales para secante y bisección
    valor_inferior = 0.02
    valor_superior = 0.022


    g_inferior = funcion_numerica(valor_inferior)
    g_superior = funcion_numerica(valor_superior)

    # justificacion de los valores
    print("\nJustificación de los valores iniciales")
    print("-------------------------------------")

    print("\nSecante:")
    print("Se seleccionan x0 = 0.02 y x1 = 0.022 porque en la gráfica")
    print("se observa que la raíz positiva se encuentra entre valores")
    print("cercanos a este intervalo.")
    print("El método de la secante utiliza dos aproximaciones iniciales")
    print("y no requiere calcular la derivada de la función.")

    print("\nBisección:")
    print("Se selecciona el intervalo [0.02, 0.022] porque se observa")
    print("un cambio de signo de g(f) dentro de este intervalo.")
    print("g(0.02) =", g_inferior)
    print("g(0.022) =", g_superior)
    print("g(0.02)*g(0.022) =", g_inferior*g_superior)

    if g_inferior*g_superior < 0:
        print("Como g(0.02)*g(0.022) < 0, el intervalo cumple")
        print("la condición necesaria para aplicar Bisección.")
    else:
        print("El intervalo no cumple la condición necesaria para Bisección.")


if __name__ == "__main__":

    # Datos de la tubería.
    diametro = 0.25
    rugosidad = 0.00015
    reynolds = 120000

    # Función de Colebrook escrita como g(f) = 0
    g_str = "1/sqrt(x) + 2*log(0.00015/(3.7*0.25) + 2.51/(120000*sqrt(x)))/log(10)"

    analisis_funcion(g_str)

    iter_max = 1000
    tol = 1e-8
    resultados = []


    # Secante

    start_time = time.perf_counter()
    xk, erk, k, conv = secante(g_str, 0.02, 0.022, iter_max, tol)
    total_time = time.perf_counter() - start_time

    resultados.append(["Secante", xk, erk, k, total_time, conv])


    # Bisección

    start_time = time.perf_counter()
    xk, erk, k, conv = biseccion(g_str, 0.02, 0.022, iter_max, tol)
    total_time = time.perf_counter() - start_time

    resultados.append(["Bisección", xk, erk, k, total_time, conv])


    # Newton-Raphson



    # Steffensen
 


    # Falsa Posición
   


    # Müller
   


    print("\nResultados parciales")
    print("--------------------")

    print(f"\n{'Método':<12} | {'xk':<14} | {'erk':<12} | {'k':<5} | {'Tiempo (ms)':<12} | {'conv':<5}")
    print("-" * 72)

    for fila in resultados:
        metodo, xk, erk, k, t_seg, conv = fila
        tiempo_ms = t_seg*1000
        print(f"{metodo:<12} | {xk:<14.10f} | {erk:<12.4e} | {k:<5} | {tiempo_ms:<12.4f} | {conv:<5}")


    # Gráfica comparativa de errores

    metodos = []
    errores = []

    for fila in resultados:
        metodo, xk, erk, k, t_seg, conv = fila
        metodos.append(metodo)
        errores.append(erk)

    plt.figure(figsize=(10, 6))
    plt.bar(metodos, errores)

    plt.yscale("log")
    plt.title("Comparación del error final de los métodos")
    plt.xlabel("Método")
    plt.ylabel("Error final")
    plt.grid(True)
    plt.show()
    