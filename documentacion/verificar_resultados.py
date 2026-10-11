# Comprobar los resultados guardados sin repetir el costoso calculo de Jacobi.
from pathlib import Path
import importlib.util
import json
import sys

import numpy as np
from PIL import Image


def verificar_resultados():
    # Importar las funciones reales de la Parte III desde la carpeta del proyecto.
    documentacion = Path(__file__).resolve().parent
    programa = documentacion.parent / "Portafolio - Grupo 10"
    sys.path.insert(0, str(programa))
    especificacion = importlib.util.spec_from_file_location("parte3", programa / "parte3.py")
    parte3 = importlib.util.module_from_spec(especificacion)
    especificacion.loader.exec_module(parte3)
    datos = programa / "dataset_att" / "dataset_att"
    resultados = json.loads((documentacion / "resultados.json").read_text(encoding="utf-8"))

    # Reconstruir el entrenamiento para verificar su preparacion independientemente.
    S, etiquetas, rutas = parte3.cargar_entrenamiento(datos)
    promedio = S.mean(axis=1)
    A = S - promedio[:, None]
    G = A @ A.T
    assert S.shape == (2576, 360) and G.shape == (2576, 2576)
    assert np.array_equal(etiquetas, np.repeat(np.arange(1, 41), 9))
    assert 0 <= S.min() <= S.max() <= 1
    for indice, ruta in enumerate(rutas):
        with Image.open(ruta) as imagen:
            esperado = np.array(imagen.convert("L"), dtype=float).reshape(-1, order="F") / 255
        assert np.array_equal(S[:, indice], esperado)
    assert np.max(np.abs(A.mean(axis=1))) < 1e-12
    assert np.array_equal(G, G.T)

    # Revisar todos los pares propios, no solo el indicador de convergencia.
    with np.load(documentacion / "pares_propios.npz") as pares:
        valores, vectores = pares["valores"], pares["vectores"]
    assert valores.shape == (2576,) and vectores.shape == (2576, 2576)
    assert np.all(np.isfinite(valores)) and np.all(np.isfinite(vectores))
    residuos = np.linalg.norm(G @ vectores - vectores * valores, axis=0)
    ortogonalidad = np.linalg.norm(vectores.T @ vectores - np.eye(2576), "fro")
    assert residuos.max() / np.linalg.norm(G, "fro") < 1e-7
    assert ortogonalidad < 1e-7
    Uk, r, k, variacion = parte3.seleccionar_caras_base(valores, vectores)
    assert r == resultados["r"] <= 359 and k == resultados["k"]
    retenidos = np.sort(valores[valores > 1e-10 * valores.max()])[::-1]
    anterior = retenidos[:k-1].sum() / retenidos.sum()
    assert anterior < 0.95 <= variacion

    # Confirmar que las coordenadas y el modelo guardado corresponden al entrenamiento.
    with np.load(documentacion / "modelo.npz") as modelo:
        assert np.array_equal(modelo["promedio"], promedio)
        assert np.array_equal(modelo["Uk"], Uk)
        assert np.allclose(modelo["X"], Uk.T @ A, rtol=1e-12, atol=1e-12)
        assert np.array_equal(modelo["etiquetas"], etiquetas)
        X = modelo["X"].copy()
    originales = [elemento.copy() for elemento in (promedio, Uk, X)]

    # Repetir las cuarenta consultas con el modelo fijo y comprobar cada vecino.
    assert len(resultados["predicciones"]) == 40
    for persona, fila in enumerate(resultados["predicciones"], start=1):
        consulta = datos / f"s{persona}" / "10.png"
        assert consulta not in rutas and consulta.is_file()
        predicha, cercana, distancia = parte3.identificar_fotografia(
            consulta, promedio, Uk, X, etiquetas, rutas)
        assert fila["real"] == persona and fila["predicha"] == predicha
        assert fila["archivo"] == f"s{persona}/10.png"
        assert fila["cercana"] == cercana.relative_to(datos).as_posix()
        assert fila["correcta"] == (persona == predicha)
        assert np.isclose(fila["distancia"], distancia, rtol=1e-12, atol=1e-12)
    assert all(np.array_equal(actual, original)
               for actual, original in zip((promedio, Uk, X), originales))
    aciertos = sum(fila["correcta"] for fila in resultados["predicciones"])
    assert resultados["aciertos"] == aciertos
    assert resultados["exactitud"] == 100 * aciertos / 40

    # Comprobar el desempate exacto con una proyeccion de diagnostico independiente.
    base_prueba = np.eye(2576, 1)
    coordenadas_iguales = np.zeros((1, 360))
    _, cercana, _ = parte3.identificar_fotografia(
        datos / "s1" / "10.png", promedio, base_prueba,
        coordenadas_iguales, etiquetas, rutas)
    assert cercana == rutas[0]
    comprobaciones = {"verificaciones": "aprobadas", "residuo_maximo": float(residuos.max()),
                      "ortogonalidad": float(ortogonalidad), "r": r, "k": k,
                      "variacion": float(variacion), "variacion_anterior": float(anterior),
                      "aciertos": aciertos, "exactitud": 100 * aciertos / 40}
    (documentacion / "verificaciones.json").write_text(
        json.dumps(comprobaciones, indent=2), encoding="utf-8")
    print(comprobaciones)


# Ejecutar comprobaciones solo cuando este archivo se invoca directamente.
if __name__ == "__main__":
    verificar_resultados()
