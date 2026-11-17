"""Regenera todos los resultados LOCALES verificables, sin exigir Internet.

A = comparación de convenciones + importación opcional de CSV OpenCOR.
B = metodología de procesado probada con DEMO SINTÉTICA;
para adquirir registros verdaderos usar 'src/allen_lote.py' en equipo conectado.
"""
import json
import platform
import subprocess
import sys
from pathlib import Path
import numpy,scipy,matplotlib

ROOT=Path(__file__).resolve().parent
LOG=ROOT/'resultados'/'validacion'
LOG.mkdir(parents=True,exist_ok=True)
commands=[
 ('hh', [sys.executable,'src/hh_fourier.py']),
 ('corriente', [sys.executable,'src/estudio_parametrico.py']),
 ('gibbs', [sys.executable,'src/comparacion_gibbs.py']),
 ('equivalencia', [sys.executable,'src/validacion_physiome.py']),
 ('sensibilidad', [sys.executable,'src/experimentos_adicionales.py']),
 ('demostracion_no_real', [sys.executable,'src/analisis_experimental.py','--demo-sintetica']),
 ('pruebas', [sys.executable,'-m','unittest','discover','-s','tests','-v'])]
res=[]
for label,cmd in commands:
    print('Ejecutando',label,flush=True)
    p=subprocess.run(cmd,cwd=ROOT,capture_output=True,text=True)
    (LOG/f'{label}_stdout.txt').write_text(p.stdout,encoding='utf8')
    (LOG/f'{label}_stderr.txt').write_text(p.stderr,encoding='utf8')
    res.append({'etapa':label,'comando':' '.join(cmd),'codigo_salida':p.returncode,
                'estado':'OK' if p.returncode==0 else 'ERROR'})
    if p.returncode!=0:
        print(p.stderr[:1000],file=sys.stderr)
        break
meta={'python':sys.version,'platform':platform.platform(),
      'numpy':numpy.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,
      'pasos':res,'allen_descargado_en_este_pipeline':False,
      'physiome_CellML_ejecutado_en_este_pipeline':False}
(LOG/'trazabilidad_ejecucion.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False),encoding='utf8')
print(json.dumps(res,indent=2,ensure_ascii=False))
if any(x['codigo_salida']!=0 for x in res):sys.exit(1)
