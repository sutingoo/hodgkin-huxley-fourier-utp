# Hodgkin–Huxley + Fourier | UTP · Series y Transformadas

Proyecto académico reproducible de simulación biofísica, descomposición y reconstrucción por series/FFT, análisis de armónicos y evaluación con registros electrofisiológicos del Allen Institute. Comparación numérica independiente con un modelo CellML ejecutado en OpenCOR.

> **Estado:** Se dispone de código de HHlab, 10 pruebas automáticas (última ejecución documentada), resultados sintéticos, 1.400 evaluaciones armónicas de 7 ventanas de 2 especímenes de ratón, auditoría local de 2 NWB informada por el autor y comparación OpenCOR–Python para un mismo protocolo. NO se incluyen los archivos originales NWB ni las series temporales experimentales por barrido; deben descargarse desde AllenSDK. El CellML utilizado es el del tutorial de OpenCOR; no debe atribuirse automáticamente al experimento de Physiome e/672.

## Inicio rápido (Ubuntu 24.04)
```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python ejecutar_todo.py
python -m unittest discover -s tests -v
python comparar_opencor.py
```
`comparar_opencor.py` usa las fuentes `data/physiome/` y guarda resultados en `resultados/physiome_externo/`. La ejecución completa no requiere conexión para los experimentos sintéticos y OpenCOR exportado.

## Allen Institute — adquirir los registros originales
Se necesita Python 3.10 aparte (AllenSDK 2.16.2 puede fallar con Python 3.12):
```bash
python -m pip install uv
uv python install 3.10
uv venv --seed --python 3.10 .venv-allen
source .venv-allen/bin/activate
python -m pip install -r requirements-allen.txt
python src/allen_lote.py --specimens 324257146 320654829 --max-sweeps 5
deactivate
source .venv/bin/activate
python verificar_allen.py
python src/analisis_experimental.py
```
La disponibilidad de los barridos en la API y la descarga están sujetas a la conexión y a las versiones de sus dependencias. `verificar_allen.py` compara hashes respecto al manifiesto local y no sustituye una verificación por parte de Allen.

## Qué archivos consultar
- `src/hh_fourier.py`: EDO, integración, FFT y métricas.
- `src/allen_lote.py`: obtención del conjunto experimental.
- `src/analisis_experimental.py`: reconstrucción de cada ventana experimental.
- `comparar_opencor.py`: comparación independiente OpenCOR ↔ HHlab.
- `resultados/resumen_validacion.json`: métricas del caso matemático.
- `resultados/experimental/metricas_por_ventana.csv`: 1400 filas / 7 barridos / 2 métodos / N=1..100.
- `resultados/experimental/figuras/`: gráficas Allen.
- `resultados/physiome_externo/`: tablas y gráficas de comparación de OpenCOR.
- `docs/TRAZABILIDAD.md`: procedencia, limitaciones y fuentes.
- `docs/GITHUB_PUBLICACION.md`: cómo publicar y citar un commit en Word.

## Interpretación científica breve
- Hodgkin–Huxley a corriente constante de 10 µA/cm²: periodo ≈14.636 ms; 25 armónicos: RMSE≈0.450 mV, energía AC retenida≈99.9645 %. Primer N que cumple los criterios académicos: 21.
- Registros Allen: con N=25, el RMSE medio FFT convencional es 0.3436 mV (5 ventanas espécimen 320654829), 1.9316 mV (2 ventanas espécimen 324257146).
- El puente lineal reduce RMSE en N=25 en las siete ventanas; su coste de descripción adicional debe incluirse en una comparación real de compresión.
- OpenCOR tutorial vs HHlab, mismo pulso homologado: RMSE en Run #1≈0.0002193 mV. Las ejecuciones repetidas NO son observaciones biológicas.

## Validez y licencias
Código académico y conjuntos de datos con procedencias distintas. No afirmar ahorro de energía física, utilidad clínica ni óptimo universal. Antes de añadir licencia al repositorio, acordar la titularidad con el equipo. Los datos originales NWB provienen del Allen Institute; cita siempre su documentación y el identificador de espécimen.

## Informe académico editable
`docs/Informe_Final_UTP_HH.docx` es el informe integrado para el curso. **Antes de la entrega**, sustituir su marca `PENDIENTE DE PUBLICACIÓN` por la URL real de este repositorio y añadir el commit SHA usado.
