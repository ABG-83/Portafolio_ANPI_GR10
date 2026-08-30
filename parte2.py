from parte1 import secante, biseccion, newton_raphson, steffensen, falsa_posicion, muller

import sympy as sp
import numpy as np
import time
import matplotlib.pyplot as plt

def analisis_funcion(f_str):
    """
    Para justificar los valores iniciales
    """

    x = sp.symbols("x")
    funcion_simbolica = sp.sympify(f_str)
    funcion_numerica = sp.lambdify(x, funcion_simbolica, "numpy")
    
    
    x_vals = np.linspace(-3, 3, 1000)
    y_vals = funcion_numerica(x_vals)
    
    plt.figure(figsize=(10, 6))
    plt.plot(x_vals, y_vals, 'b-', linewidth=2, label='g(x)')
    plt.axhline(y=0, color='r', linestyle='--', alpha=0.5, label='g(x) = 0')
    plt.axvline(x=0, color='k', linestyle='-', alpha=0.3)
    plt.grid(True, alpha=0.3)
    plt.xlabel('x', fontsize=12)
    plt.ylabel('g(x)', fontsize=12)
    plt.title('Análisis de la función g(x) = ln(x²+1) - cos(x) - 4x² + 10', fontsize=14)
    plt.legend()
    plt.show()
    
    # Valores utilizados para Secante y Bisección
    valor_1 = funcion_numerica(1)
    valor_2 = funcion_numerica(2)
    
    
    print("""
    - Secante:
        Se seleccionan x0=1 y x1=2 porque en la gráfica se observa una
        raíz positiva entre estos valores. El método de la secante requiere
        dos aproximaciones iniciales y no necesita calcular la derivada.
        
    - Bisección:
        Se selecciona el intervalo [1, 2] porque en la gráfica se observa
        una raíz positiva dentro de este intervalo.
          """)
    
    print(" f(1) =", valor_1)
    print(" f(2) =", valor_2)
    print(" f(1)*f(2) =", valor_1*valor_2)
    
    if valor_1*valor_2 < 0:
        print(" Como f(1)*f(2) < 0, existe cambio de signo en [1, 2].")
        print(" El intervalo cumple la condición necesaria para Bisección.")
    else:
        print(" El intervalo no cumple la condición necesaria para Bisección.")

    
    print("""
    - Falsa posición:
        Mirando la gráfica se nota que la raíz positiva está entre [1, 2]
        Además se cumple el teorema de Bolzano
        
    - Müller:
        Se toman x1=1 y x2=2, valores que encierran la raíz
        Se usa x0=0, valor que permite generar correctamente la parábola 
        trazada por los 3 puntos
          """)
    
if __name__ == "__main__":
    
    f_str = "log(x**2 + 1) - cos(x) - 4*x**2 + 10"
    analisis_funcion(f_str)
    
    iter_max = 1000
    tol = 1e-8
    
    resultados = []
    
    # secante
    start_time = time.perf_counter()
    xk, erk, k, conv = secante(f_str, 1, 2, iter_max, tol)
    total_time = time.perf_counter() - start_time
    
    resultados.append(["Secante", xk, erk, k, total_time, conv])
    
    # biseccion
    
    start_time = time.perf_counter()
    xk, erk, k, conv = biseccion(f_str, 1, 2, iter_max, tol)
    total_time = time.perf_counter() - start_time
    
    resultados.append(["Bisección", xk, erk, k, total_time, conv])
    
    # newton raphson
    
    start_time = time.perf_counter()
    xk, erk, k, conv = newton_raphson(f_str, 1, iter_max, tol)
    total_time = time.perf_counter() - start_time
    
    resultados.append(["Newton-Raphson", xk, erk, k, total_time, conv])
    
    # steffensen
    
    start_time = time.perf_counter()
    xk, erk, k, conv = steffensen(f_str, 1, iter_max, tol)
    total_time = time.perf_counter() - start_time
    
    resultados.append(["Steffensen", xk, erk, k, total_time, conv])
    
    # falsa posicion
    start_time = time.perf_counter()
    xk, erk, k, conv = falsa_posicion(f_str, 1, 2, iter_max, tol)
    total_time = time.perf_counter() - start_time
    
    resultados.append(["Falsa Posición", xk, erk, k, total_time, conv])
    
    # muller
    start_time = time.perf_counter()
    xk, erk, k, conv = muller(f_str, 0, 1, 2, iter_max, tol)
    total_time = time.perf_counter() - start_time
    
    resultados.append(["Müller", xk, erk, k, total_time, conv])
    


    print(f"\n{'Método':<16} | {'xk':<12} | {'erk':<12} | {'k':<5} | {'Tiempo (ms)':<12} | {'conv':<5}")
    print("-" * 75)
    for fila in resultados:
        metodo, xk, erk, k, t_seg, conv = fila
        tiempo_ms = t_seg * 1000
        print(f"{metodo:<16} | {xk:<12.8f} | {erk:<12.4e} | {k:<5} | {tiempo_ms:<12.4f} | {conv:<5}")
       
        
       
    # Gráfica comparativa de errores obtenidos
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
    
    
    
    # Gráfica comparativa de tiempos 
    metodos = []
    tiempos_ms = []
    
    for fila in resultados:
        metodo, xk, erk, k, t_seg, conv = fila
        metodos.append(metodo)
        tiempos_ms.append(t_seg*1000)
    
    plt.figure(figsize=(10, 6))
    plt.bar(metodos, tiempos_ms)
    
    plt.title("Comparación del tiempo de ejecución de los métodos")
    plt.xlabel("Método")
    plt.ylabel("Tiempo de ejecución (ms)")
    plt.grid(True)
    plt.show()
    
    
    
    
    
    
    
    
    
    