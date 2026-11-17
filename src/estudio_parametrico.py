"""Evalúa si el N recomendado cambia al variar I_ext (datos sintéticos)."""
import csv
import sys
from pathlib import Path
import numpy as np
from hh_fourier import simulate, isolated_cycle, fft_coefficients, reconstruct, metrics

ROOT = Path(__file__).resolve().parents[1]
recs=[]
for current in [0.0, 6.0, 8.0, 10.0, 12.0, 15.0]:
    sol = simulate(current=current, duration_ms=250)
    try:
        t, y, T, peaks = isolated_cycle(sol)
        c = fft_coefficients(y[0]); v=y[0]
        n_opt=None
        for n in range(1,101):
            m=metrics(t,v,reconstruct(c,n),c,n,T)
            if m['RMSE_mV']<1 and m['energia_AC_retenida_pct']>=99.9 and m['error_pico_mV']<1:
                n_opt=n;break
        recs.append({'I_ext_uA_cm2':current,'estado':'oscilación periódica',
                    'periodo_ms':T,'frecuencia_Hz':1000/T,
                    'N_min_criterio':n_opt,'RMSE_N25_mV':metrics(t,v,reconstruct(c,25),c,25,T)['RMSE_mV']})
    except ValueError:
        recs.append({'I_ext_uA_cm2':current,'estado':'no se observaron 3 espigas periódicas',
                     'periodo_ms':'','frecuencia_Hz':'','N_min_criterio':'','RMSE_N25_mV':''})
p=ROOT/'resultados'/'estudio_parametrico.csv'
with p.open('w',newline='',encoding='utf8') as f:
    writer=csv.DictWriter(f,fieldnames=list(recs[0]));writer.writeheader();writer.writerows(recs)
for rec in recs:print(rec)
