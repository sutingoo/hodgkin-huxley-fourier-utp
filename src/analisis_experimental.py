"""Análisis Fourier de clips de electrofisiología (Allen o demo SINTÉTICA).

Sin fuga de procedencia: solo se etiqueta como REAL si el CSV está en el
manifiesto de una descarga de AllenSDK. Los resultados de demo no se mezclan.
No se estima frecuencia de disparo neuronal a partir de un clip de 12 ms:
1/T aquí es solo separación de bins FFT de la ventana artificial.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from hh_fourier import fft_coefficients, reconstruct, ac_energy_fraction

ROOT=Path(__file__).resolve().parents[1]


def read_trace(path):
    tab=np.genfromtxt(path,delimiter=',',names=True,encoding='utf8')
    if tab.ndim==0:raise ValueError('CSV sin muestras suficientes')
    cols=tab.dtype.names
    if 'tiempo_ms' not in cols or 'voltaje_mV' not in cols:
        raise ValueError('Encabezados CSV esperados: tiempo_ms, voltaje_mV')
    t=np.asarray(tab['tiempo_ms'],dtype=float)
    v=np.asarray(tab['voltaje_mV'],dtype=float)
    keep=np.isfinite(t)&np.isfinite(v)
    t,v=t[keep],v[keep]
    if len(v)<100 or np.any(np.diff(t)<=0):raise ValueError('Ventana insuficiente o tiempo inválido')
    dt=np.diff(t)
    if (np.max(dt)-np.min(dt))/np.median(dt)>0.01:raise ValueError('FFT exige muestreo uniforme; remuestrear con anti-aliasing si procede')
    # irfft en hh_fourier requiere longitud par: se recorta el punto final.
    if len(v)%2: v=v[:-1];t=t[:-1]
    return t,v


def single_trace(path, original_tipo='NO_VERIFICADO',tag='desconocido'):
    t,v=read_trace(path)
    M=len(v); dt=float(np.median(np.diff(t)))
    T=M*dt  # la ventana dura M pasos, NO t[-1]-t[0]
    pp=float(np.ptp(v))
    # Corrige discontinuidad en los bordes de la extensión periódica.
    trend=np.linspace(v[0],v[-1],M)
    modes={'FFT_cruda':(v,None),'FFT_residuo_con_puente':(v-trend,trend)}
    all_rows=[]
    examples={'t':t,'v':v}
    for mode,(wave,bridge) in modes.items():
        c=fft_coefficients(wave)
        efrac_den=float(np.var(wave))
        for n in range(1,min(100,M//2-1)+1):
            r=reconstruct(c,n)
            if bridge is not None:r=r+bridge
            diff=r-v
            energy=100*ac_energy_fraction(c,n) if efrac_den>1e-18 else 100.
            row={'id':tag,'tipo':original_tipo,'metodo':mode,'N':n,
                 'duracion_ventana_ms':T,'muestras':M,'fs_Hz':1000/dt,
                 'separacion_bines_Hz':1000/T,'salto_extremos_mV':float(v[-1]-v[0]),
                 'RMSE_mV':float(np.sqrt(np.mean(diff*diff))),
                 'NRMSE_pct_Vpp':float(np.sqrt(np.mean(diff*diff))/pp*100),
                 'error_pico_mV':float(abs(np.max(r)-np.max(v))),
                 'error_max_mV':float(np.max(np.abs(diff))),
                 'energia_residuo_AC_retenida_pct':energy,
                 'criterio_ok':int(np.sqrt(np.mean(diff*diff))<1 and abs(np.max(r)-np.max(v))<1
                                  and energy>=99.9)}
            all_rows.append(row)
            if n in (5,10,20,25,50):examples[mode+'_N'+str(n)]=r
    return all_rows,examples


def analyze(input_dir, out, allow_unverified=False):
    out.mkdir(parents=True,exist_ok=True)
    provenance=input_dir/'procedencia_allen.json'
    real_set=set()
    if provenance.is_file():
        obj=json.loads(provenance.read_text(encoding='utf8'))
        for row in obj.get('registros',[]):
            if row.get('tipo')=='EXPERIMENTAL_REAL':real_set.add(Path(row['csv']).name)
    csvs=sorted(input_dir.glob('allen_*_sweep_*.csv'))
    if allow_unverified:csvs=sorted(set(csvs)|set(input_dir.glob('*.csv')))
    if not csvs:raise SystemExit('No hay archivos CSV Allen: primero ejecutar allen_lote.py con conexión, o prueba demo explícita.')
    results=[];failures=[]
    for path in csvs:
        tipo='ALLEN_REAL_DECLARADO_EN_MANIFIESTO' if path.name in real_set else 'NO_VERIFICADO'
        if 'DEMO' in path.name.upper() or 'SINTETIC' in path.name.upper():tipo='SINTETICO_DEMOSTRACION'
        if tipo=='NO_VERIFICADO' and not allow_unverified:
            failures.append({'archivo':path.name,'motivo':'No consta en procedencia_allen.json'});continue
        try:
            rows,examples=single_trace(path,tipo,path.stem)
            results+=rows
            t,v=examples['t'],examples['v']
            fig,ax=plt.subplots(figsize=(10,4.5))
            ax.plot(t,v,'k-',lw=1.7,label='Original')
            for name in ['FFT_cruda_N25','FFT_residuo_con_puente_N25']:
                if name in examples:ax.plot(t,examples[name],lw=1.1,label=name)
            ax.set(title=f'{path.stem} — {tipo}',xlabel='Tiempo (ms)',ylabel='Voltaje intracelular (mV)')
            ax.grid(alpha=.25);ax.legend(fontsize=8);fig.tight_layout()
            fig.savefig(out/(path.stem+'.png'),dpi=150);plt.close(fig)
        except Exception as ex:failures.append({'archivo':path.name,'motivo':str(ex)})
    if not results:raise SystemExit(f'No se pudo analizar ninguna señal. Fallos: {failures}')
    with (out/'metricas_por_ventana.csv').open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
    # NO reportar distribuciones de cortes como si los clips fueran sujetos independientes.
    groups=sorted(set(r['id'] for r in results))
    aggregate=[]
    for tag in groups:
        for mode in ['FFT_cruda','FFT_residuo_con_puente']:
            sub=[x for x in results if x['id']==tag and x['metodo']==mode]
            feasible=[x['N'] for x in sub if x['criterio_ok']]
            n25=next((x for x in sub if x['N']==25),None)
            aggregate.append({'id':tag,'tipo':sub[0]['tipo'],'metodo':mode,
                             'N_min_criterios':min(feasible) if feasible else 'NO_ALCANZADO_EN_RANGO',
                             'rmse_N25_mV':n25['RMSE_mV'] if n25 else '',
                             'peak_error_N25_mV':n25['error_pico_mV'] if n25 else '',
                             'endpoint_jump_mV':sub[0]['salto_extremos_mV']})
    with (out/'resumen_por_ventana.csv').open('w',newline='',encoding='utf8') as f:
        w=csv.DictWriter(f,fieldnames=list(aggregate[0]));w.writeheader();w.writerows(aggregate)
    report={'n_ventanas_analizadas':len(groups),'n_ventanas_declaradas_reales':sum('ALLEN_REAL' in x['tipo'] for x in aggregate)//2,
            'n_ventanas_sinteticas':sum('SINTETICO' in x['tipo'] for x in aggregate)//2,
            'n_ventanas_sin_verificar':sum(x['tipo']=='NO_VERIFICADO' for x in aggregate)//2,
            'fallos':failures,'metodologia':'12ms ventana nativa, FFT cruda y puente lineal + FFT + restauración. Mismo N para ambas.',
            'limitaciones':'Los coeficientes no equivalen a periodo de disparo. Criterios arbitrados para estudio; requiere resultados por células independientes.'}
    (out/'reporte_analisis.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf8')
    print(json.dumps(report,indent=2,ensure_ascii=False))
    return report


def demo(out):
    """Demo totalmente SINTÉTICA para verificar pipeline. NO es dato Allen."""
    out.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(20261116)
    t=np.arange(-4.,8.,0.02)
    for i,(sigma,shift) in enumerate([(0.2,0.),(.5,.2)]):
        # Pulso compuesto artificial NO generador electrofisiológico validado.
        v=-65 + 96*np.exp(-((t-shift)/.48)**2)-13*np.exp(-((t-shift-1.05)/1.15)**2)
        v += rng.normal(0,sigma,len(t))
        with (out/f'DEMO_SINTETICA_{i+1}.csv').open('w',newline='',encoding='utf8') as f:
            w=csv.writer(f);w.writerow(['tiempo_ms','voltaje_mV','estimulo_pA'])
            w.writerows(zip(t,v,np.full(len(t),120)))
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--input-dir',type=Path,default=ROOT/'data'/'experimental_allen')
    p.add_argument('--out',type=Path,default=ROOT/'resultados'/'experimental')
    p.add_argument('--demo-sintetica',action='store_true')
    p.add_argument('--allow-unverified',action='store_true')
    a=p.parse_args()
    if a.demo_sintetica:
        a.input_dir=demo(ROOT/'data'/'demo_sintetica')
        a.out=ROOT/'resultados'/'demo_sintetica'
        a.allow_unverified=True
    analyze(a.input_dir,a.out,a.allow_unverified)
