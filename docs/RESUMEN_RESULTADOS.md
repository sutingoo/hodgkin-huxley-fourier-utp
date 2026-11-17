# Resumen de resultados y alcance

## Caso periódico Hodgkin–Huxley
T = 14.63621 ms, frecuencia = 68.3237 Hz, Vmax = +30.4303 mV, Vmin = −74.8962 mV, RMSE N25 = 0.44997 mV, energía **AC** retenida N25 = 99.96454 %. N mínimo según los 3 umbrales académicos = 21.

## Registro intracelular Allen (ratón)
7 ventanas de 12 ms y 2400 muestras analizadas a 200 kHz, de 2 especímenes: 5 ventanas de 320654829 y 2 de 324257146. Para N25: media RMSE convencional 0.3436 y 1.9316 mV por espécimen; para puente 0.2577 y 1.8998 mV. Criterio académico original cumplido en 5/7 ventanas, en ambos métodos. **El puente contabiliza una tendencia adicional y su energía relativa no comparte denominador con la FFT cruda**.

## Validación OpenCOR/CellML
Archivo del tutorial, pulso +20 µA/cm² en Python equivalente a −20 µA/cm² entre 10 y 10.5 ms en CellML; V_Python = −65−V_CellML. Run #1: RMSE 0.000219327 mV sobre 4501 tiempos; máximo absoluto 0.00301741 mV. No confundir este protocolo con corriente constante de +10 µA/cm².

## Estado de publicación
Repositorio sugerido `sutingoo/hodgkin-huxley-fourier-utp` **PENDIENTE DE PUBLICACIÓN**. Aportar URL real y hash commit al informe después del push.
