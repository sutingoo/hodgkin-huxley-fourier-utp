"""Análisis Fourier en una ventana de voltaje medida en CSV.

AVISO: la FFT en ventana finita no convierte la actividad neuronal en periódica.
Es importante reportar ventana, extremos, sampling y características del experimento.

Uso: python src/analizar_csv_real.py data/experimental_allen/allen_....csv
Archivo esperado: 'tiempo_ms','voltaje_mV' (generado por obtener_datos_allen.py).
"""
import argparse
import csv
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from hh_fourier import fft_coefficients,reconstruct,metrics

ap=argparse.ArgumentParser()
ap.add_argument('csv',type=Path)
ap.add_argument('--out',type=Path,default=Path(__file__).resolve().parents[1]/'resultados'/'experimental')
args=ap.parse_args()
raw=np.genfromtxt(args.csv,delimiter=',',names=True,encoding='utf8')
t=np.asarray(raw['tiempo_ms'],dtype=float);v=np.asarray(raw['voltaje_mV'],dtype=float)
keep=np.isfinite(t)&np.isfinite(v)
t=t[keep];v=v[keep]
if len(t)<100 or np.any(np.diff(t)<=0):raise SystemExit('CSV insuficiente o tiempo no creciente')
T=float(t[-1]-t[0]); t2=np.arange(4096)*T/4096
v2=np.interp(t2,t-t[0],v)
c=fft_coefficients(v2)
args.out.mkdir(parents=True,exist_ok=True)
rows=[]
for n in [1,3,5,10,15,20,25,30,40,50,75,100]:
    rows.append(metrics(t2,v2,reconstruct(c,n),c,n,T))
with (args.out/'metricas_experimentales.csv').open('w',newline='',encoding='utf8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
fig,ax=plt.subplots(figsize=(10,4))
ax.plot(t2,v2,color='black',lw=2,label='Medición experimental')
for n in [5,10,25,50]:ax.plot(t2,reconstruct(c,n),label=f'N={n}',lw=1)
ax.set(xlabel='Tiempo relativo a la ventana (ms)',ylabel='Voltaje intracelular (mV)',
       title='Análisis espectral de registro experimental (ventana finita)')
ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(args.out/'reconstruccion_experimental.png',dpi=150)
(args.out/'ADVERTENCIA.txt').write_text(
    'Datos experimentales de otra especie y/o tipo neuronal, NO validan por sí solos HH axón de calamar.\n'
    f'Archivo fuente: {args.csv}\nVentana procesada: {T:.6f} ms.\n'
    f'Diferencia entre extremos: {v2[-1]-v2[0]:.6f} mV; si es grande habrá fuga espectral por extensión periódica.\n'
    'Los N armónicos dependen de duración de ventana; no comparar directamente con un periodo HH.\n',encoding='utf8')
print('Resultados experimentales guardados en',args.out)
