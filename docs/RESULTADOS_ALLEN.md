# Análisis científico de las curvas de Fourier en señales del Allen Institute

**Proyecto:** Series y Transformadas — Hodgkin–Huxley y reconstrucción armónica  
**Dataset analizado:** `metricas_por_ventana.csv` (archivo de salida de HHlab)  
**Naturaleza del estudio:** análisis secundario de métricas numéricas proporcionadas por el usuario, no una nueva descarga ni una nueva ejecución de AllenSDK.  
**Estado de procedencia:** el CSV declara `ALLEN_REAL_DECLARADO_EN_MANIFIESTO`. Esta etiqueta no sustituye la auditoría del NWB original, el manifiesto ni sus metadatos.

## 1. Datos y metodología

El CSV contiene **1.400 filas**, **7 ventanas** de **2 especímenes**, dos procedimientos de reconstrucción y un barrido entero de **N = 1 a 100** armónicos por ventana y método. Todas las ventanas se registran como segmentos de **12 ms** con **2.400 muestras** y **frecuencia de muestreo de 200 kHz**, por lo que la separación de bins es aproximadamente 83,333 Hz. Los 1.400 puntos de evaluación no son 1.400 observaciones biológicas independientes.

Métodos comparados:

- `FFT_cruda`: reconstrucción espectral directa de cada ventana.
- `FFT_residuo_con_puente`: se resta una tendencia lineal entre los extremos, se aplica FFT al residuo y se restaura la tendencia para evaluar el error frente a la señal original.

Para comparar las reconstrucciones se utiliza el error cuadrático medio en milivoltios:

\[
\mathrm{RMSE}(N)=\sqrt{\frac{1}{M}\sum_{i=0}^{M-1}\left[V_i-\widehat V_{N,i}\right]^2}.
\]

El CSV también informa el error de pico `abs(max(V_estimado)-max(V_original))` y el error máximo local `max(abs(V_estimado-V_original))`.

La regla operacional implementada en el proyecto es:

\[
\mathrm{criterio\_ok} = [\mathrm{RMSE}<1\,\mathrm{mV}]\land[\mathrm{error\ de\ pico}<1\,\mathrm{mV}]\land[E_{\mathrm{AC}}\ge99.9\%].
\]

**Comprobación:** las 1.400 filas concuerdan exactamente con esta regla. Los umbrales son decisiones metodológicas del proyecto, no umbrales clínicos.

> Precaución sobre energía: `energia_residuo_AC_retenida_pct` usa como denominador la energía AC de la **señal original** para FFT cruda, pero la energía AC del **residuo tras retirar el puente** para el segundo método. Los porcentajes no deben compararse directamente entre los métodos como si midieran exactamente el mismo referente. El RMSE sí está calculado frente al voltaje original en ambos métodos.

## 2. Resultados resumidos a N=25

| Espécimen | Ventanas | RMSE medio FFT cruda (mV) | RMSE medio con puente (mV) | Reducción relativa de medias |
|:--|--:|--:|--:|--:|
| 320654829 | 5 | 0,3436 | 0,2577 | 25,01 % |
| 324257146 | 2 | 1,9316 | 1,8998 | 1,65 % |

Los RMSE de los **siete** barridos se reducen levemente o considerablemente con el puente lineal en N=25. El método convencional cumple la regla para 5/7 ventanas, al igual que el puente: son las cinco del espécimen 320654829. En el espécimen 324257146, las dos ventanas fallan con N=25 para ambos métodos.

| Espécimen / barrido | N mínimo FFT | N mínimo puente | RMSE FFT N25 (mV) | RMSE puente N25 (mV) | Error pico FFT N25 (mV) |
|:--|--:|--:|--:|--:|--:|
| 320654829 / 51 | 20 | 18 | 0,4290 | 0,2637 | 0,5356 |
| 320654829 / 52 | 20 | 18 | 0,4108 | 0,2645 | 0,5221 |
| 320654829 / 53 | 19 | 18 | 0,3357 | 0,2638 | 0,5079 |
| 320654829 / 54 | 19 | 18 | 0,3066 | 0,2609 | 0,5184 |
| 320654829 / 56 | 17 | 17 | 0,2361 | 0,2355 | 0,2939 |
| 324257146 / 33 | 61 | 61 | 1,7324 | 1,7156 | 4,6867 |
| 324257146 / 35 | 59 | 61 | 2,1308 | 2,0840 | 4,8862 |

**Interpretación:** N=25 no es suficiente para las dos ventanas de un espécimen bajo los criterios actuales. Tampoco existe evidencia para asumir que los cortes anteriores se generalizan a neuronas humanas, de ratón o a diferentes configuraciones de electrodos sin revisar metadatos y ampliar el muestreo por célula.

## 3. Cuestiones detectadas en la auditoría

### 3.1. El criterio de error máximo no está incluido

Para N=25, el error máximo absoluto es mayor de 1 mV en las 7 ventanas y en ambos métodos. Por ejemplo, el barrido 51 del espécimen 320654829 cumple los tres criterios originales con FFT cruda, pero tiene un error máximo local de **3,8313 mV**.

Evaluamos adicionalmente el requisito exploratorio \(\max |V-\hat V|<1\,\mathrm{mV}\), sin reemplazar los criterios del código. En N=1..100:

- FFT cruda solo alcanza los cuatro criterios para una de siete ventanas (barrido 56, N=31).
- FFT con puente los alcanza en siete de siete ventanas, con primer N entre 31 y 91.

Este resultado describe una **prueba de sensibilidad de criterios**, no prueba una ventaja clínica ni de compresión. Una validación más exigente debería analizar también ancho de pulso, errores temporales y detección de eventos.

### 3.2. El primer N que cumple no siempre significa cumplimiento sostenido

En la ventana 324257146/sweep 35 con FFT cruda, el primer N válido es **59**, pero **N=60 vuelve a fallar** por error de pico (1,027 mV), y N=61 vuelve a pasar (0,956 mV). Por tanto, el criterio no es estrictamente monótono. La columna `N_min_criterios` representa el **primer N aceptable**, no necesariamente un umbral que garantice aceptación para todo N mayor.

### 3.3. Número de armónicos y frecuencia de muestreo

Para ventanas de 12 ms, \(\Delta f = 1/T \approx 83,333\,\mathrm{Hz}\), de modo que N=25 incluye componentes hasta ~2.083,33 Hz y N=61 hasta ~5.083,33 Hz. Estas frecuencias corresponden a **bins de una ventana**, no a la frecuencia de disparo neuronal. La señal fue originalmente adquirida a 200 kHz, y los resultados no demuestran una reducción directa equivalente en la frecuencia de muestreo ni en consumo energético.

### 3.4. Independencia estadística, alcance biológico y validación externa

Se han analizado 5 barridos de un espécimen y 2 del otro. Las repeticiones por célula están correlacionadas; **n biológico = 2**, no 7. No se calculan p-valores, intervalos de confianza por célula ni inferencia poblacional con este archivo. El archivo de métricas no incluye especies, tipos celulares, archivo NWB, SHA-256 o condiciones completas de estímulo; se requiere acompañarlo con `procedencia_allen.json` y, cuando sea posible, el NWB original.

Este análisis no sustituye la comparación directa del modelo HH de Python contra el CellML oficial de Physiome/OpenCOR: es otro componente del protocolo científico.

## 4. Conclusiones sustentables

1. En las siete ventanas, el puente lineal produce un RMSE menor que la FFT cruda a N=25.
2. Las cinco ventanas del espécimen 320654829 cumplen los criterios originales con N=25; las dos del espécimen 324257146 no cumplen.
3. El primer N aceptable oscila entre 17 y 61, dependiendo de la ventana y el método; el resultado no es universal.
4. El error máximo local revela una limitación significativa del criterio de aceptación vigente.
5. Aún es necesario ampliar a múltiples **células independientes**, verificar procedencia y metadatos y ejecutar la validación de Physiome.

## 5. Productos incluidos

- `metricas_por_ventana.csv`: CSV original adjunto por el usuario, sin alteraciones.
- `resumen_ventanas.csv`: 7 filas con métricas emparejadas, cortes y errores locales.
- `resumen_especimenes.csv`: agregación **descriptiva** por espécimen.
- `curvas_por_especimen.csv`: valores promedio por N, espécimen y método (no inferencia poblacional).
- `auditoria_numerica.json`: controles automáticos y trazabilidad de los resúmenes.
- `figuras/`: gráficos de convergencia RMSE, comparación N=25, cortes mínimos, error local y no monotonía del error de pico.
- `analizar_csv.py`: script para regenerar estos productos sin AllenSDK (requiere NumPy y Matplotlib).

## 6. Fuente de datos y documentación técnica

**Archivo primario de esta entrega:** `metricas_por_ventana.csv` adjunto en esta conversación, etiquetado por HHlab como `ALLEN_REAL_DECLARADO_EN_MANIFIESTO`. La verificación independiente del origen exige los NWB y el manifiesto que no forman parte de este archivo.

- AllenSDK, documentación de datos electrofisiológicos: https://alleninstitute.github.io/AllenSDK/cell_types.html
- Allen Brain Map: https://celltypes.brain-map.org/
- Modelo Hodgkin–Huxley, Physiome Model Repository: https://models.cellml.org/e/672
- Hodgkin & Huxley (1952), *The Journal of Physiology*: https://doi.org/10.1113/jphysiol.1952.sp004764
