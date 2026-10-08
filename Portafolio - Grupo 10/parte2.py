import time
import numpy as np
import matplotlib.pyplot as plt

from parte1 import metodo_qr
# from parte1 import metodo_jacobi


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
    print("ANÁLISIS COMPARATIVO")
    print("=" * 80)

    for resultado in resultados:
        print(f"\nMétodo {resultado['nombre']}")
        print(f"  Error final:              {resultado['erk']:.6e}")
        print(f"  Iteraciones:              {resultado['k']}")
        print(f"  Tiempo:                   {resultado['tiempo']:.6f} s")
        print(f"  Convergencia:             {resultado['conv']}")
        print(f"  Residuo máximo:            {resultado['r_max']:.6e}")
        print(f"  Error de ortogonalidad:    {resultado['er_ort']:.6e}")

    if len(resultados) >= 2:
        primero = resultados[0]
        segundo = resultados[1]

        diferencia = np.linalg.norm(
            primero["lambdas"] - segundo["lambdas"],
            ord=2
        )

        print("\nComparación de valores propios:")
        print(
            f"  Norma de la diferencia entre "
            f"los valores propios: {diferencia:.6e}"
        )

        print("\nNota:")
        print(
            "  Los signos de los vectores propios pueden cambiar "
            "entre métodos sin que esto implique una diferencia "
            "en el vector propio."
        )


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
        # ("Jacobi", metodo_jacobi),
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