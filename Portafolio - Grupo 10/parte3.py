
# Importar las herramientas de datos, imagenes y el metodo existente.
from pathlib import Path
import time

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from parte1 import metodo_qr, metodo_jacobi
from parte2 import construir_matriz, evaluar_metodo, mostrar_tabla


#  Carga del entrenamiento 
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


# Seleccionar componentes exclusivamente a partir del espectro de entrenamiento.
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


def comprobar_pares_propios(G, valores, vectores, error, iteraciones, conv, tiempo):
    # Separar el criterio de parada del metodo de la calidad de sus pares propios.
    if valores.shape != (2576,) or vectores.shape != (2576, 2576):
        raise ValueError("Dimensiones incorrectas en los pares propios de G.")
    if not np.all(np.isfinite(valores)) or not np.all(np.isfinite(vectores)):
        raise ValueError("El metodo devolvio pares propios no finitos.")
    residuos = np.linalg.norm(G @ vectores - vectores * valores, axis=0)
    residuo_maximo = float(residuos.max())
    residuo_relativo = residuo_maximo / np.linalg.norm(G, "fro")
    ortogonalidad = np.linalg.norm(vectores.T @ vectores - np.eye(2576), "fro")
    error_normas = np.max(np.abs(np.linalg.norm(vectores, axis=0) - 1))
    print(f"\nJacobi: error de parada = {error:.6e}; barridos = {iteraciones}; conv = {conv}")
    print(f"Tiempo del metodo: {tiempo:.2f} s ({tiempo/60:.2f} min)")
    print(f"Residuo maximo de los 2576 pares: {residuo_maximo:.6e}")
    print(f"Residuo maximo / norma de Frobenius de G: {residuo_relativo:.6e}")
    print(f"Error de ortogonalidad: {ortogonalidad:.6e}")
    print(f"Error maximo de norma de los vectores: {error_normas:.6e}")

    # Estos limites son controles locales adicionales, no la tolerancia de parada.
    if not conv or not np.isfinite(error) or error >= 1e-5:
        raise RuntimeError("Jacobi no alcanzo el criterio de parada exigido.")
    if residuo_relativo > 1e-7 or ortogonalidad > 1e-7 or error_normas > 1e-7:
        raise RuntimeError("Los pares propios no superaron los controles de calidad.")
    if valores.max() <= 0 or valores.min() < -1e-10 * valores.max():
        raise RuntimeError("El espectro calculado contiene valores incompatibles con G.")


def mostrar_caras_base(promedio, valores, vectores, carpeta_figuras):
    # Mostrar las primeras seis caras retenidas aunque k pudiera ser menor que seis.
    orden = np.argsort(valores)[::-1]
    conservar = valores[orden] > 1e-10 * valores.max()
    caras = vectores[:, orden][:, conservar]
    figura, ejes = plt.subplots(2, 4, figsize=(8, 6))
    ejes = ejes.ravel()
    ejes[0].imshow(promedio.reshape(56, 46, order="F"), cmap="gray", vmin=0, vmax=1)
    ejes[0].set_title("Rostro promedio")
    for indice in range(min(6, caras.shape[1])):
        ejes[indice+1].imshow(caras[:, indice].reshape(56, 46, order="F"), cmap="gray")
        ejes[indice+1].set_title(f"Cara base {indice+1}")
    for eje in ejes:
        eje.axis("off")
    figura.tight_layout()
    figura.savefig(carpeta_figuras / "caras_base.png", dpi=200)


def evaluar_pruebas(carpeta_datos, promedio, Uk, X, etiquetas, rutas, carpeta_figuras):
    # Evitar modificaciones accidentales del modelo durante las cuarenta consultas.
    for matriz in (promedio, Uk, X, etiquetas):
        matriz.setflags(write=False)
    predicciones = []
    print("\nArchivo       Real  Predicha  Acerto")
    for persona in range(1, 41):
        ruta = carpeta_datos / f"s{persona}" / "10.png"
        if ruta in rutas:
            raise ValueError("Una fotografia de prueba aparece en el entrenamiento.")
        predicha, cercana, distancia = identificar_fotografia(
            ruta, promedio, Uk, X, etiquetas, rutas)
        correcta = persona == predicha
        predicciones.append((ruta, persona, predicha, correcta, cercana, distancia))
        archivo = ruta.relative_to(carpeta_datos).as_posix()
        print(f"{archivo:<13} {persona:>4} {predicha:>9} {'Si' if correcta else 'No':>7}")

    # Informar hechos medidos sin atribuir los errores a una causa no comprobada.
    aciertos = sum(fila[3] for fila in predicciones)
    print(f"\nAciertos: {aciertos}/40; exactitud: {100*aciertos/40:.2f}%")
    confusiones = [fila for fila in predicciones if not fila[3]]
    print(f"Confusiones reales: {len(confusiones)}")
    for ruta, persona, predicha, _, cercana, distancia in confusiones:
        print(f"{ruta.relative_to(carpeta_datos).as_posix()}: persona {persona} -> "
              f"persona {predicha}; vecino {cercana.relative_to(carpeta_datos).as_posix()}; "
              f"distancia {distancia:.6f}.")
    if confusiones:
        print("Hecho: en cada confusion otra identidad quedo mas cerca en el espacio reducido.")
        print("Las parejas permiten observar diferencias de expresion, postura e intensidades.")
        print("Iluminacion, expresion y perdida de informacion son posibles causas; "
              "estas pruebas no demuestran cual produjo cada error.")
    print("Conservar el 95 % de variacion no garantiza un 95 % de exactitud.")
    print("Se identifica entre 40 personas conocidas; no se evalua rechazo de desconocidos.")

    # Presentar las cuarenta parejas en cuatro figuras con etiquetas legibles.
    for pagina in range(4):
        figura, ejes = plt.subplots(5, 4, figsize=(8, 10))
        for posicion, fila in enumerate(predicciones[pagina*10:(pagina+1)*10]):
            ruta, persona, predicha, correcta, cercana, _ = fila
            renglon, columna = posicion // 2, (posicion % 2)*2
            for desplazamiento, (imagen_ruta, titulo) in enumerate(
                    [(ruta, "Prueba"), (cercana, "Vecino")]):
                eje = ejes[renglon, columna+desplazamiento]
                with Image.open(imagen_ruta) as imagen:
                    eje.imshow(np.array(imagen.convert("L")), cmap="gray", vmin=0, vmax=255)
                archivo = imagen_ruta.relative_to(carpeta_datos).as_posix()
                eje.set_title(f"{titulo}\n{archivo}", fontsize=11,
                              color="black" if correcta else "red")
                eje.axis("off")
        figura.tight_layout()
        figura.savefig(carpeta_figuras / f"parejas_{pagina+1}.png", dpi=180)


def parte3():
    # Comparar la version actual de ambos metodos sin depender de resultados guardados.
    matriz_comparacion = construir_matriz(200)
    resultados_comparacion = []
    for nombre, metodo in [("QR", metodo_qr), ("Jacobi", metodo_jacobi)]:
        resultado = evaluar_metodo(nombre, metodo, matriz_comparacion, 10000, 1e-8)
        resultados_comparacion.append(resultado)
    mostrar_tabla(resultados_comparacion)
    qr, jacobi = resultados_comparacion
    print("\nSeleccion: Jacobi, con la comparacion real de esta ejecucion.")
    print(f"Residuos QR/Jacobi: {qr['r_max']:.6e} / {jacobi['r_max']:.6e}; "
          f"ortogonalidad: {qr['er_ort']:.6e} / {jacobi['er_ort']:.6e}.")
    print("La tabla compara tambien error de parada, iteraciones y tiempo; "
          "conv por si solo no garantiza pares propios validos.")
    # Si cambia la evidencia que sustenta Jacobi, exigir revisar la seleccion.
    if (not jacobi["conv"] or jacobi["r_max"] > qr["r_max"]
            or jacobi["er_ort"] > 1e-7 or jacobi["r_max"] > 1e-7):
        raise RuntimeError("La comparacion actual no respalda Jacobi; revise la eleccion.")

    # Resolver las rutas junto al script para poder ejecutarlo desde cualquier carpeta.
    carpeta_programa = Path(__file__).resolve().parent
    carpeta_datos = carpeta_programa / "dataset_att" / "dataset_att"
    if not carpeta_datos.is_dir():
        raise FileNotFoundError(f"No se encontro el dataset: {carpeta_datos}")
    carpeta_figuras = carpeta_programa / "figuras"
    carpeta_figuras.mkdir(exist_ok=True)

    # Construir el promedio y el centrado exclusivamente con las 360 imagenes conocidas.
    S, etiquetas, rutas = cargar_entrenamiento(carpeta_datos)
    promedio = np.mean(S, axis=1)
    A = S - promedio[:, None]
    G = A @ A.T

    # Comprobar preparacion y correspondencia antes del calculo numerico costoso.
    assert S.shape == (2576, 360) and A.shape == (2576, 360)
    assert G.shape == (2576, 2576) and promedio.shape == (2576,)
    assert np.all(np.isfinite(S)) and 0 <= S.min() <= S.max() <= 1
    assert np.array_equal(etiquetas, np.repeat(np.arange(1, 41), 9))
    assert all(ruta == carpeta_datos / f"s{persona}" / f"{foto}.png"
               for ruta, (persona, foto) in zip(
                   rutas, [(persona, foto) for persona in range(1, 41) for foto in range(1, 10)]))
    assert all((carpeta_datos / f"s{persona}" / "10.png").is_file() for persona in range(1, 41))
    with Image.open(rutas[0]) as imagen:
        original = np.array(imagen.convert("L"), dtype=float) / 255
    reconstruida = S[:, 0].reshape(56, 46, order="F")
    error_reconstruccion = np.max(np.abs(original - reconstruida))
    error_centrado = np.max(np.abs(A.mean(axis=1)))
    assert error_reconstruccion == 0 and error_centrado < 1e-12
    print("\nEntrenamiento: 360 fotografias; prueba: 40 fotografias reservadas.")
    print("Dimensiones S, promedio, A y G:", S.shape, promedio.shape, A.shape, G.shape)
    print(f"Intensidades normalizadas: [{S.min():.8f}, {S.max():.8f}]")
    print(f"Reconstruccion: {error_reconstruccion:.6e}; centrado: {error_centrado:.6e}")
    print("Orden verificado: s1/1.png a s40/9.png; 10.png excluida del entrenamiento.")

    # Visualizar el preprocesamiento conservando la escala de la imagen centrada.
    figura, ejes = plt.subplots(1, 3, figsize=(8, 4))
    imagenes = [reconstruida, promedio.reshape(56, 46, order="F"),
                A[:, 0].reshape(56, 46, order="F")]
    for indice, titulo in enumerate(["Original: s1/1.png", "Rostro promedio", "Imagen centrada"]):
        if indice < 2:
            ejes[indice].imshow(imagenes[indice], cmap="gray", vmin=0, vmax=1)
        else:
            ejes[indice].imshow(imagenes[indice], cmap="gray")
        ejes[indice].set_title(titulo)
        ejes[indice].axis("off")
    figura.tight_layout()
    figura.savefig(carpeta_figuras / "preprocesamiento.png", dpi=160)

    # Medir el costo y solicitar el avance interno del metodo existente.
    inicio = time.perf_counter()
    print("\nCalculando Jacobi sobre G completa; puede tardar decenas de minutos.", flush=True)
    print("El porcentaje y el tiempo restante corresponden al barrido actual, no al total.")
    print("La estimacion puede variar; el numero de barridos hasta converger es desconocido.")
    valores, vectores, error, iteraciones, conv = metodo_jacobi(
        G, 10**6, 1e-5, mostrar_progreso=True)
    tiempo = time.perf_counter() - inicio
    comprobar_pares_propios(G, valores, vectores, error, iteraciones, conv, tiempo)

    # Seleccionar el minimo k solo con entrenamiento, antes de abrir las consultas.
    Uk, r, k, variacion = seleccionar_caras_base(valores, vectores)
    retenidos = np.sort(valores[valores > 1e-10 * valores.max()])[::-1]
    variacion_anterior = np.sum(retenidos[:k-1]) / np.sum(retenidos)
    if r > 359 or not (variacion_anterior < 0.95 <= variacion):
        raise RuntimeError("La seleccion de las caras base no supero las comprobaciones.")
    X = proyectar_entrenamiento(A, Uk)
    assert Uk.shape == (2576, k) and X.shape == (k, 360)
    print(f"\nr = {r}; k = {k}; variacion conservada = {variacion:.8%}")
    print(f"Variacion con k-1 = {variacion_anterior:.8%}")
    print("Dimensiones Uk y X:", Uk.shape, X.shape)
    mostrar_caras_base(promedio, valores, vectores, carpeta_figuras)
    evaluar_pruebas(carpeta_datos, promedio, Uk, X, etiquetas, rutas, carpeta_figuras)

    # Guardar solo figuras; mostrarlas al final evita interrumpir el entrenamiento.
    print("\nComprobaciones y evaluacion completadas. Figuras guardadas en:", carpeta_figuras)
    plt.show()
    plt.close("all")


# Permitir importar las funciones sin ejecutar automaticamente el entrenamiento.
if __name__ == "__main__":
    parte3()


