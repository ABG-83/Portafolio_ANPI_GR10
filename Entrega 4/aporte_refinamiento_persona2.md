# Refinamiento de malla: aporte de Persona 2

Se ejecutaron los tres esquemas con tf = 10 s.

## Refinamiento espacial

| Esquema | dx (m) | dt (s) | mu | E_f (C) | Tiempo (s) | Residuo maximo |
| --- | --- | --- | --- | --- | --- | --- |
| Euler explicito | 0.040000 | 0.002500 | 0.015625 | 3.503491e-02 | 0.014470 | no aplica |
| Euler explicito | 0.020000 | 0.002500 | 0.062500 | 6.051987e-03 | 0.021403 | no aplica |
| Euler explicito | 0.010000 | 0.002500 | 0.250000 | 1.210226e-03 | 0.050224 | no aplica |
| Euler implicito | 0.040000 | 0.002500 | 0.015625 | 4.227195e-02 | 1.023050 | 2.778176e-16 |
| Euler implicito | 0.020000 | 0.002500 | 0.062500 | 1.331059e-02 | 1.987402 | 1.841674e-16 |
| Euler implicito | 0.010000 | 0.002500 | 0.250000 | 6.050195e-03 | 3.961223 | 2.215171e-16 |
| Crank-Nicolson | 0.040000 | 0.002500 | 0.015625 | 3.865365e-02 | 1.159042 | 1.845314e-16 |
| Crank-Nicolson | 0.020000 | 0.002500 | 0.062500 | 9.681517e-03 | 2.084391 | 1.894668e-16 |
| Crank-Nicolson | 0.010000 | 0.002500 | 0.250000 | 2.420211e-03 | 4.058983 | 1.911460e-16 |

Euler explicito: los errores de las tres mallas, en el orden de la tabla, son 3.503491e-02, 6.051987e-03 y 1.210226e-03 C; los cocientes E_anterior/E_siguiente son 5.789 y 5.001. El error disminuye en los dos refinamientos. Se reduce dx a la mitad manteniendo dt = 0.0025 s. El error espacial es de segundo orden; un cociente cercano a cuatro es compatible con su predominio, pero el error temporal fijo y su posible compensacion impiden exigir ese cociente al error total.

Euler implicito: los errores de las tres mallas, en el orden de la tabla, son 4.227195e-02, 1.331059e-02 y 6.050195e-03 C; los cocientes E_anterior/E_siguiente son 3.176 y 2.200. El error disminuye en los dos refinamientos. Se reduce dx a la mitad manteniendo dt = 0.0025 s. El error espacial es de segundo orden; un cociente cercano a cuatro es compatible con su predominio, pero el error temporal fijo y su posible compensacion impiden exigir ese cociente al error total.

Crank-Nicolson: los errores de las tres mallas, en el orden de la tabla, son 3.865365e-02, 9.681517e-03 y 2.420211e-03 C; los cocientes E_anterior/E_siguiente son 3.993 y 4.000. El error disminuye en los dos refinamientos. Se reduce dx a la mitad manteniendo dt = 0.0025 s. El error espacial es de segundo orden; un cociente cercano a cuatro es compatible con su predominio, pero el error temporal fijo y su posible compensacion impiden exigir ese cociente al error total.

## Refinamiento temporal

| Esquema | dx (m) | dt (s) | mu | E_f (C) | Tiempo (s) | Residuo maximo |
| --- | --- | --- | --- | --- | --- | --- |
| Euler explicito | 0.020000 | 0.016000 | 0.400000 | 1.356129e-02 | 0.003856 | no aplica |
| Euler explicito | 0.020000 | 0.008000 | 0.200000 | 1.935956e-03 | 0.006986 | no aplica |
| Euler explicito | 0.020000 | 0.004000 | 0.100000 | 3.873818e-03 | 0.013361 | no aplica |
| Euler implicito | 0.020000 | 0.016000 | 0.400000 | 3.289382e-02 | 0.302245 | 2.608148e-16 |
| Euler implicito | 0.020000 | 0.008000 | 0.200000 | 2.129159e-02 | 0.606403 | 2.265746e-16 |
| Euler implicito | 0.020000 | 0.004000 | 0.100000 | 1.548759e-02 | 1.276783 | 1.924842e-16 |
| Crank-Nicolson | 0.020000 | 0.016000 | 0.400000 | 9.675555e-03 | 0.350961 | 2.204233e-16 |
| Crank-Nicolson | 0.020000 | 0.008000 | 0.200000 | 9.680138e-03 | 0.628473 | 1.822953e-16 |
| Crank-Nicolson | 0.020000 | 0.004000 | 0.100000 | 9.681284e-03 | 1.245048 | 2.053615e-16 |

Euler explicito: los errores de las tres mallas, en el orden de la tabla, son 1.356129e-02, 1.935956e-03 y 3.873818e-03 C; los cocientes E_anterior/E_siguiente son 7.005 y 0.500. El error no disminuye en ambos refinamientos. Con M = 50 fijo, el error espacial permanece. Para el modo senoidal, el operador espacial discreto amortigua menos que el continuo y Euler explicito introduce mayor amortiguamiento temporal. Ambos errores pueden compensarse; al reducir dt esa compensacion cambia, lo que explica una posible subida del error aun cumpliendo mu <= 0.5.

Euler implicito: los errores de las tres mallas, en el orden de la tabla, son 3.289382e-02, 2.129159e-02 y 1.548759e-02 C; los cocientes E_anterior/E_siguiente son 1.545 y 1.375. El error disminuye en los dos refinamientos. Con M = 50 fijo, Euler implicito tiene error temporal de primer orden. Para este dato senoidal, tanto el error espacial como el temporal reducen el enfriamiento respecto de la solucion exacta. Al reducir dt disminuye la contribucion temporal, pero permanece el error espacial; por eso los cocientes del error total no tienen que ser dos.

Crank-Nicolson: los errores de las tres mallas, en el orden de la tabla, son 9.675555e-03, 9.680138e-03 y 9.681284e-03 C; los cocientes E_anterior/E_siguiente son 1.000 y 1.000. El error no disminuye en ambos refinamientos. Con M = 50 fijo, Crank-Nicolson tiene error temporal de segundo orden. Un cociente cercano a cuatro solo se espera si domina ese error; si los errores cambian poco al reducir dt, el error espacial fijo limita la mejora. La compensacion de contribuciones tambien puede alterar los cocientes.

Euler implicito: todos los sistemas de las seis mallas cumplieron el criterio; el mayor residuo normalizado fue 2.778176e-16, inferior a 1e-10.

Crank-Nicolson: todos los sistemas de las seis mallas cumplieron el criterio; el mayor residuo normalizado fue 2.204233e-16, inferior a 1e-10.

Los tiempos de las tablas incluyen la construccion de matrices y el avance numerico; excluyen errores, tablas y graficas. Euler implicito reutiliza la factorizacion de Thomas y no tiene iteraciones internas. Al refinar el tiempo aumentan los pasos (625, 1250 y 2500); al refinar el espacio aumentan las incognitas interiores (24, 49 y 99). Las mediciones pueden fluctuar y no constituyen por si solas una estimacion del orden de complejidad.
