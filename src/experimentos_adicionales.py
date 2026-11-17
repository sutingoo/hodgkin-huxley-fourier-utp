"""Dos experimentos suplementarios reproducibles: muestreo y robustez a ruido.
Solo usa la señal Hodgkin-Huxley simulada. No es evidencia de hardware.
"""
import csv
import json
from pathlib import Path
import numpy as np
from scipy.signal import resample
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from hh_fourier import simulate, isolated_cycle, fft_coefficients, reconstruct

ROOT=Path(__file__).resolve().parents[1]

def save_rows(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with open(path,'w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def run():
    out=ROOT/'resultados'/'sensibilidad';out.mkdir(parents=True,exist_ok=True)
    t,y,T,_=isolated_cycle(simulate(duration_ms=220))
    v=y[0];rng=np.random.default_rng(20261116)
    noiserows=[]
    for sigma in [0,0.5,1.0,2.0]:
        for rep in range(40):
            measured=v+rng.normal(0,sigma,len(v))
            c=fft_coefficients(measured)
            for N in [5,10,15,20,21,25,30,40,50,75,100]:
                for filt in ['truncada','Lanczos']:
                    r=reconstruct(c,N,lanczos=(filt=='Lanczos'))
                    noiserows.append({'sigma_ruido_simulado_mV':sigma,'repeticion':rep,
                                      'N':N,'filtro':filt,
                                      'rmse_vs_original_limpia_mV':float(np.sqrt(np.mean((r-v)**2))),
                                      'rmse_vs_observada_ruidosa_mV':float(np.sqrt(np.mean((r-measured)**2))),
                                      'error_pico_vs_limpia_mV':float(abs(np.max(r)-np.max(v)))})
    save_rows(out/'montecarlo_ruido.csv',noiserows)
    agg=[]
    for sigma in [0,0.5,1.0,2.0]:
        for N in [5,10,15,20,21,25,30,40,50,75,100]:
            for filt in ['truncada','Lanczos']:
                a=np.array([z['rmse_vs_original_limpia_mV'] for z in noiserows if
                            z['sigma_ruido_simulado_mV']==sigma and z['N']==N and z['filtro']==filt])
                agg.append({'sigma_mV':sigma,'N':N,'filtro':filt,'promedio_RMSE_mV':float(a.mean()),
                            'p2_5_RMSE_mV':float(np.quantile(a,.025)),
                            'p97_5_RMSE_mV':float(np.quantile(a,.975))})
    save_rows(out/'resumen_montecarlo.csv',agg)
    fig,ax=plt.subplots(figsize=(10,5))
    for sigma in [.5,1.,2.]:
        rows=[z for z in agg if z['sigma_mV']==sigma and z['filtro']=='truncada']
        ax.plot([z['N'] for z in rows],[z['promedio_RMSE_mV'] for z in rows],'-o',label=f'ruido σ={sigma} mV')
    ax.set(xlabel='N armónicos',ylabel='RMSE contra referencia sin ruido (mV)',title='Compromiso armónicos-ruido (simulado, n=40 por nivel)')
    ax.legend();ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'ruido_vs_N.png',dpi=160);plt.close(fig)
    # Estudio de muestreo sin filtro antialias predefinido. Tras remuestreo se compara al ideal.
    sample_rows=[]
    for fs_requested in [2000,3500,5000,10000,20000,50000]:
        M=int(np.round(T*fs_requested/1000))
        # M equiespaciados sobre T; fs efectivo = M*1000/T
        tnew=np.arange(M)*T/M
        vnew=np.interp(tnew,t,v,period=T)
        interpolated=resample(vnew,len(v)) # extensión periódica asumida explícitamente
        sample_rows.append({'f_s_solicitado_Hz':fs_requested,'n_muestras_ciclo':M,
                            'f_s_efectivo_Hz':1000*M/T,'frecuencia_nyquist_Hz':500*M/T,
                            'rmse_recuperacion_sin_antialias_mV':float(np.sqrt(np.mean((interpolated-v)**2))),
                            'error_max_mV':float(np.max(np.abs(interpolated-v))),
                            'error_pico_mV':float(abs(np.max(interpolated)-np.max(v)))})
    save_rows(out/'muestreo_sintetico.csv',sample_rows)
    fig,ax=plt.subplots(figsize=(9,4.5))
    ax.plot([x['f_s_efectivo_Hz'] for x in sample_rows],
            [x['rmse_recuperacion_sin_antialias_mV'] for x in sample_rows],'o-')
    ax.set(xlabel='Frecuencia de muestreo efectiva (Hz)',ylabel='RMSE reconstrucción (mV)',title='Error de muestreo sobre HH sintético (sin anti-alias explícito)')
    ax.grid(alpha=.2);fig.tight_layout();fig.savefig(out/'muestreo.png',dpi=160);plt.close(fig)
    summary={'tipo':'SIMULACION_SINTETICA','seed':20261116,'repeticiones_por_sigma':40,
             'niveles_ruido_mV':[0,.5,1,2],
             'advertencia':'No se midió ruido fisiológico ni consumo; interpolación de Fourier asume señal periódica y no sustituye diseño ADC antialias.',
             'muestreo':sample_rows,
             'rmse_N25_ruido_1mV_promedio':next(x['promedio_RMSE_mV'] for x in agg if x['sigma_mV']==1 and x['N']==25 and x['filtro']=='truncada')}
    (out/'resumen_sensibilidad.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf8')
    print('Sensibilidad creada:',out,'RMSE N25 ruido1=',summary['rmse_N25_ruido_1mV_promedio'])
    return summary

if __name__=='__main__':run()
