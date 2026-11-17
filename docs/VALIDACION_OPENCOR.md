# Validación numérica cruzada del modelo Hodgkin–Huxley: Python/HHlab frente a OpenCOR/CellML

## Objetivo

Comprobar que la implementación del sistema de Hodgkin–Huxley en HHlab reproduce las soluciones obtenidas con OpenCOR cuando ambos programas utilizan el **mismo protocolo eléctrico, las mismas variables de estado iniciales, los mismos parámetros y las mismas unidades**.

## Datos y trazabilidad

Los dos archivos analizados fueron proporcionados por el usuario: `hodgkin_huxley_squid_axon_model_1952_original.cellml` y `hodgkin_huxley_squid_axon_model_1952_original_data.csv`. El CSV contiene 4.501 tiempos (0–45 ms, paso 0,01 ms) y 13 columnas de ejecución por variable de estado. Se escogió `Run #1` como evaluación primaria y se revisaron las otras doce ejecuciones como comprobaciones complementarias. El archivo del modelo incluye un encabezado de generación mediante COR en 2009; **no se ha demostrado que sea exactamente la versión modular e/672 alojada en Physiome**.

## Homologación matemática de las convenciones

Para las ecuaciones del modelo CellML adjunto, el potencial de membrana `V_CellML` emplea una referencia invertida respecto a la convención habitual de voltaje absoluto de HHlab, y un potencial de reposo cercano a cero. La transformación usada en la comparación es:

\[
V_{\mathrm{HHlab}}(t)=-65\ \mathrm{mV}-V_{\mathrm{CellML}}(t)
\]

Para igualar el signo del estímulo:

\[
I_{\mathrm{HHlab}}(t)=-i_{\mathrm{Stim,CellML}}(t)
\]

El archivo CellML define un pulso de `−20 μA/cm²` entre 10 y 10,5 ms, equivalente a `+20 μA/cm²` con signo HHlab. El estímulo fuera de ese intervalo es cero. Las condiciones iniciales utilizadas son:

\[
V(0)=-65\ \mathrm{mV},\quad m(0)=0.05,\quad h(0)=0.6,\quad n(0)=0.325.
\]

Se conserva el sistema implementado por HHlab:

\[
C_m\frac{dV}{dt}=I_\mathrm{ext}-g_\mathrm{Na}m^3h(V-E_\mathrm{Na})-g_\mathrm{K}n^4(V-E_\mathrm{K})-g_L(V-E_L)
\]

con las tres ecuaciones cinéticas clásicas `dx/dt = α_x(V)(1−x)−β_x(V)x`, para `x∈{m,h,n}`. El código fuente importa **sin reescribir** la función `hh_rhs` de HHlab y resuelve el problema mediante `scipy.integrate.solve_ivp`, método DOP853, `rtol=2×10⁻¹¹`, `atol=10⁻¹³`, paso máximo de 0,05 ms. Se integra por tramos para tratar la discontinuidad de entrada en 10 y 10,5 ms.

## Indicadores de error

Para las `M=4501` muestras:

\[
\mathrm{RMSE}=\sqrt{\frac{1}{M}\sum_{i=1}^{M}(V_\mathrm{HHlab}(t_i)-V_\mathrm{OpenCOR,convertido}(t_i))^2}
\]

\[
E_\mathrm{max}=\max_i|V_\mathrm{HHlab}(t_i)-V_\mathrm{OpenCOR,convertido}(t_i)|
\]

## Resultados (Run #1)

| Indicador | Resultado |
|---|---:|
| Número de muestras | 4.501 |
| RMSE del voltaje | 0,0002193267 mV (0,2193 µV) |
| Error absoluto máximo | 0,003017413 mV (3,0174 µV) |
| Sesgo medio | −1,12 × 10⁻⁷ mV |
| Pico máximo OpenCOR convertido | +39,499177 mV |
| Pico máximo HHlab | +39,499105 mV |
| Tiempo del máximo en ambos | 12,07 ms |
| Mínimo de OpenCOR convertido | −76,176902 mV |
| Mínimo HHlab | −76,176902 mV |
| RMSE de compuerta m | 1,903 × 10⁻⁶ |
| RMSE de compuerta h | 5,165 × 10⁻⁷ |
| RMSE de compuerta n | 5,432 × 10⁻⁷ |

Para las trece ejecuciones exportadas: RMSE máximo `0,000321585 mV`, RMSE mediano `0,000147310 mV` y error absoluto máximo global `0,004364842 mV`.

## Figuras de referencia

1. `resultados/01_comparacion_voltaje.png`: superposición OpenCOR y HHlab bajo condiciones equivalentes.
2. `resultados/02_error_temporal.png`: diferencia temporal `V_HHlab−V_OpenCOR`.
3. `resultados/03_compuertas.png`: comparación de `m(t)`, `h(t)` y `n(t)`.
4. `resultados/04_rmse_por_ejecucion.png`: variación de RMSE entre las ejecuciones exportadas.

## Interpretación científica

La correspondencia observada entre ambos solucionadores, con errores de voltaje inferiores a 0,005 mV en las ejecuciones exportadas, indica que la implementación numérica de HHlab reproduce el comportamiento dinámico previsto por el modelo CellML adjunto cuando se igualan todas las condiciones. El resultado valida la coherencia matemática y computacional de esta implementación. No demuestra validación fisiológica frente a experimentos con neuronas reales, ni prueba desempeño de compresión o ahorro de energía de hardware.

## Limitaciones y próximos experimentos

- El modelo tutorial de OpenCOR utiliza un pulso breve, mientras que el experimento periódico del informe original de la UTP utiliza `I_ext = 10 µA/cm²` constante. Los resultados deben presentarse como protocolos diferentes.
- Las 13 ejecuciones proceden del mismo modelo, no de 13 datos biológicos distintos.
- El análisis de Fourier sobre registros reales del Allen Institute debe reportarse por separado, con más células y criterios explícitos.
- Para citar Physiome e/672 específicamente, deben confirmarse el identificador y el archivo modelo exactos y ejecutarse en condiciones homologadas.

### Texto de conclusión sugerido para la investigación

> Se realizó una validación cruzada de la implementación numérica del modelo de Hodgkin–Huxley desarrollada en Python frente a los resultados exportados desde OpenCOR para un archivo CellML de referencia. Tras homologar la convención de potencial eléctrico, el estímulo y las condiciones iniciales, el RMSE de voltaje para 4501 muestras fue 2,193×10⁻⁴ mV y el error máximo absoluto fue 3,017×10⁻³ mV. La alta correspondencia obtenida respalda la correcta implementación computacional de las ecuaciones, dentro de la configuración estudiada. Esta comparación no sustituye la validación fisiológica mediante registros experimentales independientes.
