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

    # Valores utilizados para justificar Newton-Raphson y Steffensen
    derivada_simbolica = sp.diff(funcion_simbolica, x)
    derivada_numerica = sp.lambdify(x, derivada_simbolica, "numpy")
    denominador_steffensen = funcion_numerica(1 + valor_1) - valor_1
    
    
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
    - Newton-Raphson:
        Se selecciona x0=1 porque en la gráfica se observa una raíz positiva
        entre 1 y 2. Además, la derivada en x0 no es cero, por lo que la
        primera iteración del método se encuentra correctamente definida.

    - Steffensen:
        Se selecciona x0=1 porque es un valor cercano a la raíz positiva
        observada entre 1 y 2. También se verifica que el denominador
        f(x0+f(x0))-f(x0) no sea cero antes de iniciar las iteraciones.
          """)

    print(" f'(1) =", derivada_numerica(1))
    print(" Denominador inicial de Steffensen =", denominador_steffensen)

    if derivada_numerica(1) != 0:
        print(" Como f'(1) es diferente de cero, Newton-Raphson puede iniciar.")
    else:
        print(" Newton-Raphson no puede iniciar porque f'(1) es cero.")

    if denominador_steffensen != 0:
        print(" Como el denominador es diferente de cero, Steffensen puede iniciar.")
    else:
        print(" Steffensen no puede iniciar porque su denominador es cero.")

    
    print("""
    - Falsa posición:
        Mirando la gráfica se nota que la raíz positiva está entre [1, 2]
        Además se cumple el teorema de Bolzano para este intervalo como se
        demostró para bisección.
        """)
        
        
    print("""
    - Müller:
        Se usan tres aproximaciones iniciales, para las cuales se toman
        x1=1 y x2=2, valores que encierran la raíz como se mira en la gráfica
        Se usa como tercer valor x0=0, el cual permite generar correctamente 
        la parábola trazada por los tres valores iniciales.
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


    # Gráfica comparativa del número de iteraciones
    metodos = []
    iteraciones = []

    for fila in resultados:
        metodo, xk, erk, k, t_seg, conv = fila
        metodos.append(metodo)
        iteraciones.append(k)

    plt.figure(figsize=(10, 6))
    plt.bar(metodos, iteraciones)

    plt.title("Comparación del número de iteraciones de los métodos")
    plt.xlabel("Método")
    plt.ylabel("Número de iteraciones")
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


    # Análisis comparativo de los seis métodos
    resultados_convergentes = [fila for fila in resultados if fila[5] == 1]
    menor_error = min(resultados, key=lambda fila: fila[2])
    menos_iteraciones = min(resultados, key=lambda fila: fila[3])
    menor_tiempo = min(resultados, key=lambda fila: fila[4])

    aproximaciones = [fila[1] for fila in resultados_convergentes]
    diferencia_aproximaciones = max(aproximaciones) - min(aproximaciones)

    print("\nAnálisis comparativo de los resultados:")

    if len(resultados_convergentes) == len(resultados):
        print("Todos los métodos alcanzaron la tolerancia solicitada y conv = 1.")
    else:
        print("No todos los métodos alcanzaron la tolerancia solicitada.")

    print(f"Las aproximaciones convergen a la raíz positiva x = {aproximaciones[0]:.10f}.")
    print(f"La diferencia máxima entre las aproximaciones es {diferencia_aproximaciones:.4e}.")
    print(f"El menor error corresponde a {menor_error[0]}, con erk = {menor_error[2]:.4e}.")
    print(f"El menor número de iteraciones corresponde a {menos_iteraciones[0]}, con k = {menos_iteraciones[3]}.")
    print(f"El menor tiempo medido corresponde a {menor_tiempo[0]}, con {menor_tiempo[4] * 1000:.4f} ms.")

    print("""
La bisección fue el método que necesitó más iteraciones. Esto concuerda con
su convergencia lineal

La falsa posición también conserva un intervalo con cambio de signo, pero utiliza
una interpolación lineal para generar la aproximación. Por esta razón puede
requerir menos iteraciones que la bisección en este problema.

La secante alcanzó la raíz con pocas iteraciones sin calcular derivadas. Su
comportamiento fue más rápido que el de los métodos de intervalo, aunque su
convergencia depende de que los valores iniciales sean adecuados.

Newton-Raphson obtuvo una aproximación con un error muy pequeño y pocas
iteraciones. El resultado es consistente con su convergencia cuadrática cerca
de una raíz simple, aunque necesita calcular la derivada y requiere que esta no
sea cero durante las iteraciones.

Steffensen convergió hacia la misma raíz sin utilizar derivadas. Necesitó más
iteraciones que Newton-Raphson y Secante para los valores iniciales elegidos,
pero mantuvo un error inferior a la tolerancia solicitada.

Müller utilizó una interpolación cuadrática a partir de tres valores iniciales.
En este problema alcanzó la solución con muy pocas iteraciones y un error muy
pequeño, aunque cada iteración requiere operaciones más complejas que los
métodos basados en interpolación lineal.

Los tiempos de ejecución son muy pequeños y pueden cambiar ligeramente entre
ejecuciones. Por ello, deben interpretarse junto con el error y el número de
iteraciones.
""")
    
