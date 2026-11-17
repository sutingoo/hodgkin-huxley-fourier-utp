# Matriz de procedencia y límites de validación

| Resultado | Fuente | Archivo de evidencia | Limitación |
|---|---|---|---|
| HH I=10 µA/cm², T≈14.636 ms | HHlab, modelo sintético | `resultados/resumen_validacion.json` | Solo axón de calamar clásico. |
| Error HH N=25 ≈0.450 mV | HHlab FFT | `resultados/metricas_armonicos.csv` | AC y ventana de un ciclo periódica. |
| Allen, 7 ventanas y 2 especímenes | AllenSDK, registros de ratón | `data/experimental_allen/procedencia_allen.json`, `resultados/experimental/metricas_por_ventana.csv` | El repositorio NO incluye NWB ni CSV individuales; auditoría local OK informada por el autor. |
| Reglas FFT verificadas en 1400 filas | Auditoría del CSV suministrado | `resultados/experimental/auditoria_numerica.json` | La regla no incluye error máximo local; umbrales académicos. |
| OpenCOR vs HHlab RMSE≈0.0002193 mV | OpenCOR CellML tutorial, Run #1 | `data/physiome/*`, `resultados/physiome_externo/*` | No es mismo estímulo que HH I constante ni afirmación de e/672 idéntico. |

## Cadena de custodia de Allen
Manifiesto con `specimen_id`, `sweep_number`, tasa de muestreo y hash SHA-256 por NWB; los siete barridos usan dos archivos NWB. Se ha omitido la ruta absoluta local en esta copia. El verificador `verificar_allen.py` debe ejecutarse sobre los NWB y CSV descargados por el usuario.

## Qué NO se demuestra
No hay prueba en implantes, uso clínico, consumo eléctrico físico o mejora de clasificación. Solo dos especímenes independientes. No tratar siete barridos como siete neuronas. No confundir energía de la señal AC original con energía AC del residuo del puente; tampoco confundir 1/T de ventana experimental con la frecuencia neuronal de disparo.

## Fuentes primarias
- Hodgkin y Huxley (1952), DOI: 10.1113/jphysiol.1952.sp004764.
- AllenSDK: https://alleninstitute.github.io/AllenSDK/cell_types.html
- OpenCOR/CellML tutorial: https://www.cellml.org/getting-started/tutorials/opencordemo/
- Physiome e/672 (referencia teórica relacionada, no ejecución certificada del mismo archivo): https://models.cellml.org/e/672
