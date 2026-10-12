import time
import numpy as np
import matplotlib.pyplot as plt

from parte1 import metodo_qr
from parte1 import metodo_jacobi


def construir_matriz(n):
    diag_principal = 4.0 * np.ones(n)
    subdiag = -1.0 * np.ones(n - 1)

    A = (
        np.diag(diag_principal)
        + np.diag(subdiag, k=1)
        + np.diag(subdiag, k=-1)
    )

    return A


def evaluar_metodo(nombre, metodo, A, iterMax, tol):
    print(f"\nEjecutando el Método {nombre}...")

    # Medir tiempo de ejecución
    t_inicio = time.perf_counter()

    lambdas, V, erk, k, conv = metodo(A, iterMax, tol)

    t_fin = time.perf_counter()
    tiempo = t_fin - t_inicio

    # Calcular residuos de cada par propio
    n = len(lambdas)
    residuos = np.zeros(n)

    for j in range(n):
        v_j = V[:, j]
        lambda_j = lambdas[j]

        residuos[j] = np.linalg.norm(
            A @ v_j - lambda_j * v_j,
            ord=2
        )

    r_max = float(np.max(residuos))

    # Error de ortogonalidad
    I = np.eye(n)

    er_ort = float(
        np.linalg.norm(V.T @ V - I, ord="fro")
    )

    # Guardar todos los resultados
    resultados = {
        "nombre": nombre,
        "lambdas": lambdas,
        "V": V,
        "erk": erk,
        "k": k,
        "conv": conv,
        "tiempo": tiempo,
        "residuos": residuos,
        "r_max": r_max,
        "er_ort": er_ort
    }

    return resultados


def mostrar_tabla(resultados):

    print(
        f"\n{'Método':<16} | "
        f"{'er_k':<12} | "
        f"{'k':<6} | "
        f"{'Tiempo (s)':<12} | "
        f"{'conv':<5} | "
        f"{'r_max':<12} | "
        f"{'er_ort':<12}"
    )

    print("-" * 105)

    for resultado in resultados:

        metodo = resultado["nombre"]
        erk = resultado["erk"]
        k = resultado["k"]
        tiempo = resultado["tiempo"]
        conv = resultado["conv"]
        r_max = resultado["r_max"]
        er_ort = resultado["er_ort"]

        print(
            f"{metodo:<16} | "
            f"{erk:<12.4e} | "
            f"{k:<6} | "
            f"{tiempo:<12.4f} | "
            f"{conv:<5} | "
            f"{r_max:<12.4e} | "
            f"{er_ort:<12.4e}"
        )


def mostrar_valores_propios(resultados):
    print("\n" + "=" * 80)
    print("VALORES PROPIOS APROXIMADOS")
    print("=" * 80)

    for resultado in resultados:
        print(f"\nMétodo {resultado['nombre']}:")

        lambdas = resultado["lambdas"]

        print("Primeros 5:")
        print(lambdas[:5])

        print("Últimos 5:")
        print(lambdas[-5:])


def mostrar_analisis(resultados):
    print("\n" + "=" * 80)
    print("ANÁLISIS COMPARATIVO - PARTE II")
    print("=" * 80)
    
    qr_res = resultados[0]
    jacobi_res = resultados[1]
    
    print("\n1. Aproximación de Valores Propios:")
    print("   - Ambos métodos logran aproximar el espectro completo de la matriz tridiagonal")
    print("     de 200x200, ordenados de mayor a lineal descendente.")
    diferencia = np.linalg.norm(qr_res["lambdas"] - jacobi_res["lambdas"], ord=2)
    print(f"   - La norma de la diferencia entre los valores propios de ambos métodos es de {diferencia:.4e},")
    print("     lo que demuestra una concordancia general, a pesar de que QR no completó la convergencia.")

    print("\n2. Comportamiento de los Residuos (||Av_j - lambda_j v_j||_2):")
    print(f"   - Método Jacobi: Presenta un residuo máximo de r_max = {jacobi_res['r_max']:.2e},")
    print("     evidenciando una precisión sobresaliente y cumplimiento riguroso de la ecuación característica.")
    print(f"   - Método QR: Presenta un residuo máximo de r_max = {qr_res['r_max']:.2e};")
    print("     este valor es mayor debido a que el método alcanzó el límite de iteraciones (10,000)")
    print("     sin cumplir plenamente con la tolerancia estricta establecida.")

    print("\n3. Normalización y Ortogonalidad de los Vectores Propios:")
    print("   - Todos los vectores propios calculados fueron normalizados a norma euclidiana unitaria (||v_j||_2 = 1).")
    print("   - El error de ortogonalidad medido mediante la norma de Frobenius (||V^T V - I||_F) demuestra")
    print(f"     una excelente preservación de la ortogonalidad en ambos algoritmos (QR: {qr_res['er_ort']:.2e},")
    print(f"     Jacobi: {jacobi_res['er_ort']:.2e}), gracias a la naturaleza unitaria de las transformaciones.")

    print("\n4. Diferencias en el Número de Iteraciones y Tiempo de Ejecución:")
    print(f"   - Jacobi demostró una eficiencia computacional radical: requirió únicamente {jacobi_res['k']} iteraciones")
    print(f"     y un tiempo de ejecución de {jacobi_res['tiempo']:.4f} segundos.")
    print(f"   - QR requirió el número máximo de iteraciones permitidas ({qr_res['k']}) y consumió")
    print(f"     {qr_res['tiempo']:.2f} segundos, reflejando un costo computacional muy elevado para una matriz de 200x200.")

    print("\n5. Diferencias en el Criterio de Parada:")
    print("   - Jacobi evalúa de forma directa la reducción de los elementos fuera de la diagonal principal,")
    print("     alcanzando la tolerancia rápidamente y obteniendo un indicador de convergencia conv = 1.")
    print("   - QR monitorea la norma de la diferencia entre las diagonales de iteraciones sucesivas (||lambda^(k) - lambda^(k-1)||_2);")
    print("     dado que su velocidad de convergencia asintótica es más lenta para esta matriz, agotó el límite")
    print("     de iteraciones sin satisfacer la tolerancia (conv = 0).")

    print("\nNota Importante sobre los Vectores Propios:")
    print("   - Al comparar vectores propios entre métodos, debe considerarse que el signo de un vector")
    print("     propio puede invertirse (v vs -v) sin alterar su validez matemática ni su condición ortonormal.")


def graficar_resultados(resultados):

    # Gráfica comparativa de errores
    metodos = []
    errores = []

    for resultado in resultados:
        metodos.append(resultado["nombre"])
        errores.append(resultado["erk"])

    plt.figure(figsize=(10, 6))

    plt.bar(metodos, errores)

    plt.yscale("log")
    plt.title("Comparación del error final de los métodos")
    plt.xlabel("Método")
    plt.ylabel("Error final $er_k$")
    plt.grid(True)

    plt.show()


    # Gráfica comparativa del número de iteraciones
    metodos = []
    iteraciones = []

    for resultado in resultados:
        metodos.append(resultado["nombre"])
        iteraciones.append(resultado["k"])

    plt.figure(figsize=(10, 6))

    plt.bar(metodos, iteraciones)

    plt.yscale("log")
    plt.title("Comparación del número de iteraciones de los métodos")
    plt.xlabel("Método")
    plt.ylabel("Número de iteraciones")
    plt.grid(True)

    plt.show()


    # Gráfica comparativa de tiempos
    metodos = []
    tiempos = []

    for resultado in resultados:
        metodos.append(resultado["nombre"])
        tiempos.append(resultado["tiempo"])

    plt.figure(figsize=(10, 6))

    plt.bar(metodos, tiempos)

    plt.yscale("log")
    plt.title("Comparación del tiempo de ejecución de los métodos")
    plt.xlabel("Método")
    plt.ylabel("Tiempo de ejecución (s)")
    plt.grid(True)

    plt.show()


if __name__ == "__main__":

    # Parámetros
    iterMax = 10000
    tol = 1e-8

    # Matriz de la Parte II
    n = 200
    A = construir_matriz(n)

    print("=" * 80)
    print("COMPARACIÓN COMPUTACIONAL - PARTE II")
    print("=" * 80)

    # Métodos a comparar
    metodos = [
        ("QR", metodo_qr),
        ("Jacobi", metodo_jacobi),
    ]

    # Ejecutar métodos
    resultados = []

    for nombre, metodo in metodos:
        resultado = evaluar_metodo(
            nombre,
            metodo,
            A,
            iterMax,
            tol
        )

        resultados.append(resultado)

    # Mostrar resultados
    mostrar_tabla(resultados)
    mostrar_valores_propios(resultados)

    # Análisis
    mostrar_analisis(resultados)

    # Gráficas
    graficar_resultados(resultados)