# Comparar las funciones existentes sin alterar los archivos de otros integrantes.
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

# Importar la matriz y la evaluacion de la Parte II sin ejecutar su programa.
RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "Portafolio - Grupo 10"))
from parte1 import metodo_qr, metodo_jacobi
from parte2 import construir_matriz, evaluar_metodo, mostrar_tabla


def comparar():
    # Usar exactamente el problema y los parametros requeridos.
    matriz = construir_matriz(200)
    resultados = []
    for nombre, metodo in [("QR", metodo_qr), ("Jacobi", metodo_jacobi)]:
        resultado = evaluar_metodo(nombre, metodo, matriz, 10000, 1e-8)
        resultados.append(resultado)
        print(nombre, resultado["k"], resultado["r_max"], flush=True)

    # Conservar todos los pares y residuos para poder auditar la comparacion.
    carpeta = Path(__file__).resolve().parent
    mostrar_tabla(resultados)
    resumen = []
    for resultado in resultados:
        np.savez(carpeta / f"comparacion_{resultado['nombre']}.npz",
                 valores=resultado["lambdas"], vectores=resultado["V"],
                 residuos=resultado["residuos"])
        resumen.append({clave: float(resultado[clave]) for clave in
                        ["erk", "k", "conv", "tiempo", "r_max", "er_ort"]}
                       | {"nombre": resultado["nombre"]})
    registro = {"n": 200, "iterMax": 10000, "tol": 1e-8,
                "resultados": resumen,
                "parte1_sha256": hashlib.sha256(
                    (RAIZ / "Portafolio - Grupo 10/parte1.py").read_bytes()).hexdigest()}
    (carpeta / "comparacion.json").write_text(
        json.dumps(registro, indent=2), encoding="utf-8")


# Permitir reutilizar el comparador sin iniciarlo al importarlo.
if __name__ == "__main__":
    comparar()
