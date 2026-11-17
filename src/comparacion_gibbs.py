"""Distingue fenómeno de Gibbs (discontinuidad) de rizado en HH continuo."""
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from hh_fourier import fft_coefficients, reconstruct

ROOT=Path(__file__).resolve().parents[1]
raw=np.genfromtxt(ROOT/'data'/'simulados'/'hh_periodo.csv',delimiter=',',names=True)
v=raw['voltaje_mV']
square=np.where(np.sin(2*np.pi*(np.arange(len(v))+.41)/len(v))>0,1.0,-1.0)
cs=fft_coefficients(square);ch=fft_coefficients(v)
Ns=[10,25,50,100,200,400]
res=[]
for n in Ns:
    rs=reconstruct(cs,n);rh=reconstruct(ch,n)
    gibbs=max(np.max(rs)-1, -1-np.min(rs))/2*100
    serr=np.max(np.abs(rh-v))
    res.append((n,round(gibbs,5),round(float(serr),6)))
print('N,overshoot_paso_pct_salto,error_max_HH_mV')
for row in res: print(*row,sep=',')
fig,ax=plt.subplots(figsize=(11,4))
t=np.arange(len(v))/len(v)
ax.plot(t,square,color='black',label='Escalón periódico (discontinuo)')
for n in [10,25,100]:ax.plot(t,reconstruct(cs,n),label=f'Fourier N={n}',alpha=.8)
ax.set(xlabel='Tiempo normalizado (ciclos)', ylabel='Amplitud',
       title='Gibbs clásico: sobreimpulso persistente en discontinuidades')
ax.set_xlim(.42,.58);ax.set_ylim(-1.4,1.4);ax.grid(alpha=.2);ax.legend()
fig.tight_layout();fig.savefig(ROOT/'resultados'/'figuras'/'09_gibbs_escalon_ideal.png',dpi=150);plt.close(fig)

with (ROOT/'resultados'/'comparacion_gibbs.csv').open('w',encoding='utf8') as f:
    f.write('N,sobreimpulso_escalon_pct_salto,error_max_HH_mV\n')
    for row in res:f.write(','.join(map(str,row))+'\n')
