# Responsabilidades: cargar el dataset facial, preparar los datos de entrenamiento
# y presentar las comprobaciones y figuras del preprocesamiento.

# -------- Importaciones --------
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image


# -------- Carga del entrenamiento --------
def cargar_entrenamiento(carpeta_datos):
    # Mantener el orden numérico permite asociar cada imagen con su persona.
    S = np.zeros((2576, 360), dtype=float)
    etiquetas = []
    rutas = []
    columna = 0

    for persona in range(1, 41):
        for foto in range(1, 10):
            ruta = carpeta_datos / f"s{persona}" / f"{foto}.png"

            # Conservar la resolución exigida antes de preparar los píxeles.
            with Image.open(ruta) as archivo:
                imagen = np.array(archivo.convert("L"), dtype=float)

            if imagen.shape != (56, 46):
                raise ValueError(f"Dimensiones incorrectas en {ruta}")

            # Normalizar y vectorizar por columnas, como está en la guia 
            imagen = imagen / 255.0
            S[:, columna] = imagen.reshape(2576, order="F")

            # Guardar la identidad y la ruta en el mismo orden que las columnas.
            etiquetas.append(persona)
            rutas.append(ruta)
            columna += 1

    return S, np.array(etiquetas), rutas


# Ordena los valores propios y sus vectores correspondientes, elimina los valores numéricamente nulos y calcula \(r\) , despues de recibir qr y jacobi 
def seleccionar_caras_base(valores, vectores):
    # Ordenar ambos resultados juntos conserva cada par valor-vector propio.
    orden = np.argsort(valores)[::-1]
    valores = np.array(valores, dtype=float)[orden]
    vectores = np.array(vectores, dtype=float)[:, orden]

    # Evitar construir un modelo sin variación positiva.
    if valores[0] <= 0:
        raise ValueError("El mayor valor propio debe ser positivo.")

    # Descartar los valores numéricamente nulos según la guía.
    conservar = valores > 1e-10 * valores[0]
    valores = valores[conservar]
    vectores = vectores[:, conservar]
    r = len(valores)

    # Las caras base deben tener longitud uno para calcular las proyecciones.
    normas = np.linalg.norm(vectores, axis=0)
    if np.any(normas == 0):
        raise ValueError("Se recibió un vector propio de norma cero.")

    vectores = vectores / normas

    # Elegir el menor número de componentes que conserve al menos el 95 %.
    variacion_acumulada = np.cumsum(valores) / np.sum(valores)
    k = int(np.searchsorted(variacion_acumulada, 0.95)) + 1
    variacion = variacion_acumulada[k - 1]
    Uk = vectores[:, :k]

    return Uk, r, k, variacion


def proyectar_entrenamiento(A, Uk):
    # Representar cada imagen mediante sus coordenadas en las caras base.
    X = Uk.T @ A
    return X

# Identificación de una fotografía
def identificar_fotografia(ruta, promedio, Uk, X, etiquetas, rutas):
    # Preparar la consulta exactamente como las imágenes de entrenamiento.
    with Image.open(ruta) as archivo:
        imagen = np.array(archivo.convert("L"), dtype=float)

    if imagen.shape != (56, 46):
        raise ValueError(f"Dimensiones incorrectas en {ruta}")

    imagen = imagen / 255.0
    f = imagen.reshape(2576, order="F")

    # Utilizar el promedio y las caras base ya calculados con entrenamiento.
    a = f - promedio
    x = Uk.T @ a

    # Comparar la consulta con las 360 fotografías en el espacio reducido.
    diferencias = X - x[:, None]
    distancias = np.linalg.norm(diferencias, axis=0)

    # argmin devuelve la primera posición si existe un empate exacto.
    indice = int(np.argmin(distancias))
    persona_predicha = int(etiquetas[indice])
    ruta_cercana = rutas[indice]
    distancia = float(distancias[indice])

    return persona_predicha, ruta_cercana, distancia


# -------- Preprocesamiento y visualizacion --------
def parte3():
    # Buscar el dataset junto al script, independientemente de la terminal.
    carpeta_programa = Path(__file__).resolve().parent
    carpeta_datos = carpeta_programa / "dataset_att" / "dataset_att"

    if not carpeta_datos.is_dir():
        raise FileNotFoundError(
            f"No se encontro la carpeta del dataset: {carpeta_datos}"
        )

    # Construir el promedio y los datos centrados solo con entrenamiento.
    S, etiquetas, rutas = cargar_entrenamiento(carpeta_datos)
    promedio = np.mean(S, axis=1)
    A = S - promedio[:, None]

    # Preparar la matriz que recibirá el método QR o Jacobi.
    G = A @ A.T

    # Comprobar que la vectorización conserva la imagen original.
    with Image.open(rutas[0]) as archivo:
        original = np.array(archivo.convert("L"), dtype=float) / 255.0

    reconstruida = S[:, 0].reshape((56, 46), order="F")
    error_reconstruccion = np.max(np.abs(original - reconstruida))

    # El promedio de los datos centrados debe ser cercano a cero.
    error_centrado = np.max(np.abs(np.mean(A, axis=1)))

    # Presentar las dimensiones y comprobaciones de esta primera etapa.
    print("Fotografías de entrenamiento:", S.shape[1])
    print("Fotografías reservadas para prueba: 40 (no utilizadas)")
    print("Dimensiones de S:", S.shape)
    print("Dimensiones del promedio:", promedio.shape)
    print("Dimensiones de A:", A.shape)
    print("Dimensiones de G:", G.shape)
    print("Dimensiones de las etiquetas:", etiquetas.shape)

    print("Primera fotografía:", rutas[0].relative_to(carpeta_datos))
    print("Última fotografía:", rutas[-1].relative_to(carpeta_datos))

    print(f"Error de reconstrucción: {error_reconstruccion:.3e}")
    print(
        "Máximo promedio absoluto después de centrar:",
        f"{error_centrado:.3e}"
    )

    # Reconstruir las imágenes usando el mismo orden por columnas.
    imagen_promedio = promedio.reshape((56, 46), order="F")
    imagen_centrada = A[:, 0].reshape((56, 46), order="F")

    imagenes = [reconstruida, imagen_promedio, imagen_centrada]
    titulos = [
        "Original: s1/1.png",
        "Rostro promedio",
        "Imagen centrada"
    ]

    # Ajustar la escala de la diferencia para visualizar sus valores negativos.
    figura, ejes = plt.subplots(1, 3, figsize=(9, 4))

    for i in range(3):
        if i < 2:
            ejes[i].imshow(imagenes[i], cmap="gray", vmin=0, vmax=1)
        else:
            ejes[i].imshow(imagenes[i], cmap="gray")

        ejes[i].set_title(titulos[i])
        ejes[i].axis("off")

    # Guardar la figura para documentar el avance del preprocesamiento.
    figura.tight_layout()
    salida = carpeta_programa / "preprocesamiento.png"
    figura.savefig(salida, dpi=160)
    print("Figura guardada en:", salida)

    # Pendiente: importar QR o Jacobi desde parte1.py según la Parte II.
    # Aplicar el método a G con iterMax = 10**6 y tol = 10**-5.
    # Después completar las caras base, las proyecciones y la identificación.
    print("Preprocesamiento completo.")
    print("Pendiente el método de valores propios.")

    plt.show()


# Permitir importar las funciones sin ejecutar automáticamente el programa.
if __name__ == "__main__":
    parte3()
