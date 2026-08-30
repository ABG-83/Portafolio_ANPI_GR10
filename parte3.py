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

    # Valores iniciales para Newton-Raphson y Steffensen
    valor_newton = 0.02
    valor_steffensen = 0.0203
    derivada_simbolica = sp.diff(funcion_simbolica, x)
    derivada_numerica = sp.lambdify(x, derivada_simbolica, "numpy")

    g_steffensen = funcion_numerica(valor_steffensen)
    argumento_steffensen = valor_steffensen + g_steffensen
    denominador_steffensen = funcion_numerica(argumento_steffensen) - g_steffensen

    # Justificacion del valor inicial de Newton-Raphson
    print("\nNewton-Raphson:")
    print("Se selecciona x0 = 0.02 porque pertenece al intervalo fisico")
    print("[0.02, 0.022], donde la grafica y el cambio de signo indican")
    print("la presencia de una raiz. Ademas, x0 es positivo y cercano")
    print("a la solucion observada en la grafica.")
    print("g'(0.02) =", derivada_numerica(valor_newton))

    if derivada_numerica(valor_newton) != 0:
        print("Como g'(0.02) es diferente de cero, Newton-Raphson puede iniciar.")
    else:
        print("Newton-Raphson no puede iniciar porque g'(0.02) es cero.")

    # Justificacion del valor inicial de Steffensen
    print("\nSteffensen:")
    print("Se selecciona x0 = 0.0203 porque la grafica permite ubicar la")
    print("raiz cerca de este valor y se mantiene la condicion fisica f > 0.")
    print("Tambien se verifica que x0 + g(x0) sea positivo y que el")
    print("denominador de la formula de Steffensen sea diferente de cero.")
    print("x0 + g(x0) =", argumento_steffensen)
    print("g(x0 + g(x0)) - g(x0) =", denominador_steffensen)

    if argumento_steffensen > 0 and denominador_steffensen != 0:
        print("Las condiciones iniciales permiten aplicar Steffensen.")
    else:
        print("Las condiciones iniciales no permiten aplicar Steffensen.")
        
        
        
    # Justificación de intervalo para falsa posición
    print("\nFalsa posición:")
    print("Similar al caso de bisección, se selecciona el intervalo")
    print("[0.02, 0.022] porque se observa un cambio de signo de ")
    print("g(f) dentro de este intervalo.")
    print("g(0.02) =", g_inferior)
    print("g(0.022) =", g_superior)
    print("g(0.02)*g(0.022) =", g_inferior*g_superior)
    
    if g_inferior*g_superior < 0:
        print("Como g(0.02)*g(0.022) < 0, el intervalo cumple")
        print("la condición necesaria para aplicar Bisección.")
    else:
        print("El intervalo no cumple la condición necesaria para Bisección.")
    
    
    # Justificación de valores iniciales para Muller
    print("\nMüller:")
    print("El método de Müller utiliza tres aproximaciones iniciales.")
    print("Se seleccionan x1 = 0.02 y x2 = 0.022 porque en la gráfica")
    print("se observa que la raíz positiva se encuentra entre valores")
    print("cercanos a este intervalo.")
    print("Además se selecciona x0 = 0.01 ya que es un valor que permite")
    print("generar correctamente la parábola a partir de los tres valores")
    print("iniciales usados para este método.")
    


if __name__ == "__main__":

    # Datos de la tubería
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

    start_time = time.perf_counter()
    xk, erk, k, conv = newton_raphson(g_str, 0.02, iter_max, tol)
    total_time = time.perf_counter() - start_time

    resultados.append(["Newton-Raphson", xk, erk, k, total_time, conv])

    # Steffensen

    start_time = time.perf_counter()
    xk, erk, k, conv = steffensen(g_str, 0.0203, iter_max, tol)
    total_time = time.perf_counter() - start_time

    resultados.append(["Steffensen", xk, erk, k, total_time, conv])


    # Falsa Posición
    
    start_time = time.perf_counter()
    xk, erk, k, conv = falsa_posicion(g_str, 0.02, 0.022, iter_max, tol)
    total_time = time.perf_counter() - start_time

    resultados.append(["Falsa posición", xk, erk, k, total_time, conv])
   


    # Müller
    
    start_time = time.perf_counter()
    xk, erk, k, conv = muller(g_str, 0.01, 0.02, 0.022, iter_max, tol)
    total_time = time.perf_counter() - start_time

    resultados.append(["Müller", xk, erk, k, total_time, conv])
   


    # Tabla comparativa de resultados
    print("\nTabla comparativa de resultados")
    print("-------------------------------")

    print(f"\n{'Método':<16} | {'xk':<14} | {'erk':<12} | {'k':<5} | {'Tiempo (ms)':<12} | {'conv':<5}")
    print("-" * 76)

    for fila in resultados:
        metodo, xk, erk, k, t_seg, conv = fila
        tiempo_ms = t_seg*1000
        print(f"{metodo:<16} | {xk:<14.10f} | {erk:<12.4e} | {k:<5} | {tiempo_ms:<12.4f} | {conv:<5}")


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


    # Grafica comparativa de tiempos de ejecucion

    metodos = []
    tiempos_ms = []

    for fila in resultados:
        metodo, xk, erk, k, t_seg, conv = fila
        metodos.append(metodo)
        tiempos_ms.append(t_seg*1000)

    plt.figure(figsize=(10, 6))
    plt.bar(metodos, tiempos_ms)

    plt.title("Comparacion del tiempo de ejecucion de los metodos")
    plt.xlabel("Metodo")
    plt.ylabel("Tiempo de ejecucion (ms)")
    plt.grid(True)
    plt.show()


    # Analisis comparativo de los resultados

    resultados_convergentes = [fila for fila in resultados if fila[5] == 1]
    menor_error = min(resultados, key=lambda fila: fila[2])
    menos_iteraciones = min(resultados, key=lambda fila: fila[3])
    menor_tiempo = min(resultados, key=lambda fila: fila[4])

    metodos_menor_error = [fila[0] for fila in resultados if fila[2] == menor_error[2]]
    metodos_menos_iteraciones = [fila[0] for fila in resultados if fila[3] == menos_iteraciones[3]]

    print("\nAnalisis comparativo de los resultados")
    print("--------------------------------------")

    if len(resultados_convergentes) == len(resultados):
        print("Todos los metodos ejecutados alcanzaron la tolerancia y conv = 1.")
    else:
        print("No todos los metodos ejecutados alcanzaron la tolerancia solicitada.")

    if resultados_convergentes:
        aproximaciones = [fila[1] for fila in resultados_convergentes]
        diferencia_aproximaciones = max(aproximaciones) - min(aproximaciones)

        print(f"Los metodos convergentes aproximan f = {aproximaciones[0]:.10f}.")
        print(f"La diferencia maxima entre aproximaciones es {diferencia_aproximaciones:.4e}.")

    print(f"El menor error corresponde a {', '.join(metodos_menor_error)}, con erk = {menor_error[2]:.4e}.")
    print(f"La menor cantidad de iteraciones corresponde a {', '.join(metodos_menos_iteraciones)}, con k = {menos_iteraciones[3]}.")
    print(f"El menor tiempo medido corresponde a {menor_tiempo[0]}, con {menor_tiempo[4] * 1000:.4f} ms.")

    print("""
La biseccion ofrece convergencia segura porque parte de un intervalo con cambio
de signo, pero su convergencia lineal hace que normalmente requiera una mayor
cantidad de iteraciones.

La falsa posicion tambien conserva un intervalo con cambio de signo y utiliza
interpolacion lineal. Puede avanzar mas rapido que la biseccion, aunque uno de
los extremos del intervalo puede permanecer fijo durante varias iteraciones.

La secante no necesita derivadas y alcanza la solucion utilizando dos valores
iniciales. Generalmente converge mas rapido que los metodos de intervalo, pero
no conserva una garantia de convergencia por cambio de signo.

Newton-Raphson presenta convergencia cuadratica cerca de una raiz simple. En
este problema alcanza un error pequeno con pocas iteraciones, aunque necesita
calcular la derivada y requiere que esta no sea cero.

Steffensen obtiene una convergencia rapida sin calcular derivadas. En este
problema fue necesario seleccionar un valor inicial cercano a la raiz para
evitar aproximaciones negativas, ya que el modelo solo esta definido para
factores de friccion positivos.

Muller utiliza interpolacion cuadratica y tres valores iniciales. Puede alcanzar
la solucion en pocas iteraciones, aunque cada iteracion requiere operaciones
mas complejas y puede producir valores complejos en otros problemas.

Los tiempos medidos son muy pequeños y pueden variar entre ejecuciones. Por
esta razon, la eficiencia debe evaluarse junto con el error, las iteraciones,
las condiciones iniciales y la seguridad de convergencia de cada metodo.
""")
    
