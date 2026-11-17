#!/usr/bin/env python3
"""Validación EXTERNA de HHlab contra archivo exportado de OpenCOR.

Este script importa hh_rhs del proyecto del usuario: src/hh_fourier.py.
La simulación de validación emplea, ÚNICAMENTE PARA COMPARAR,
las mismas condiciones que el modelo CellML adjunto:
  V_CellML(0)=0 mV (relativo al reposo e invertido),
  m(0)=0.05, h(0)=0.6, n(0)=0.325;
  i_Stim_CellML = -20 uA/cm^2 en [10, 10.5] ms; fuera = 0.
  V_HHlab = -65 - V_CellML, I_HHlab = -i_Stim_CellML.

NO valida la ejecución de corriente constante +10 uA/cm2 de HHlab;
compara dos solucionadores del MISMO modelo y protocolo de estimulación.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
SOURCE_NAME = 'hodgkin_huxley_squid_axon_model_1952_original'
DEFAULT_CSV = BASE / 'data' / 'physiome' / (SOURCE_NAME+'_data.csv')
DEFAULT_CELLML = BASE / 'data' / 'physiome' / (SOURCE_NAME+'.cellml')


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as file:
        for block in iter(lambda: file.read(2**20), b''):
            digest.update(block)
    return digest.hexdigest()


def inspect_cellml(path: Path) -> dict:
    """Comprueba los valores exactos del modelo fuente (sin modificarlo)."""
    ns = {'c':'http://www.cellml.org/cellml/1.0#', 'm':'http://www.w3.org/1998/Math/MathML'}
    root = ET.parse(path).getroot()
    def initial(component, variable):
        for child in root.findall('c:component', ns):
            if child.attrib.get('name') == component:
                for var in child.findall('c:variable', ns):
                    if var.attrib.get('name') == variable:
                        return float(var.attrib['initial_value'])
        raise ValueError(f'Variable inicial no encontrada: {component}.{variable}')
    vals = {name:initial(comp,var) for name,comp,var in [
        ('V_rel_inicial_mV','membrane','V'),('E_R','membrane','E_R'),
        ('Cm','membrane','Cm'),('m','sodium_channel_m_gate','m'),
        ('h','sodium_channel_h_gate','h'),('n','potassium_channel_n_gate','n'),
        ('g_Na','sodium_channel','g_Na'),('g_K','potassium_channel','g_K'),
        ('g_L','leakage_current','g_L')]}
    for k,target in [('V_rel_inicial_mV',0),('E_R',0),('Cm',1),('m',.05),('h',.6),('n',.325),('g_Na',120),('g_K',36),('g_L',.3)]:
        if not np.isclose(vals[k],target,atol=1e-10):
            raise ValueError(f'Modelo CellML cambiado: {k}={vals[k]}, esperado {target}. Revisar protocolo.')
    comp = next(c for c in root.findall('c:component',ns) if c.attrib.get('name')=='membrane')
    stim_eq = None
    for eq in comp.findall('.//m:apply',ns):
        ch=list(eq)
        if len(ch)>=3 and ch[0].tag.endswith('eq') and ch[1].tag.endswith('ci') and (ch[1].text or '').strip()=='i_Stim':
            stim_eq=ch[2];break
    if stim_eq is None:
        raise ValueError('No se encontró ecuación i_Stim en CellML')
    nums=[float(el.text) for el in stim_eq.iter() if el.tag.endswith('cn')]
    if len(nums)!=4 or not np.allclose(nums,[20,10,10.5,0]):
        raise ValueError(f'Pulso del modelo CellML distinto del esperado: {nums}. No usar esta configuración.')
    vals['stimulus_cellml_uA_cm2'] = -20
    vals['stimulus_start_ms'] = 10.
    vals['stimulus_end_ms'] = 10.5
    vals['estimulo_Python_uA_cm2'] = 20
    vals['cellml_model_name'] = root.attrib.get('name')
    return vals


def load_opencor(path:Path):
    # Encabezados: membrane | V (millivolt) | Run #1 ... #13
    with path.open('r',encoding='utf-8-sig',newline='') as f:
        reader=csv.reader(f)
        header=next(reader)
        vals=np.asarray(list(reader),dtype=float)
    t=vals[:,0]
    if t.size<100 or not np.all(np.isfinite(t)) or np.any(np.diff(t)<=0):
        raise ValueError('CSV: eje temporal inválido')
    if len(header)!=vals.shape[1]:
        raise ValueError('CSV: el número de columnas no coincide')
    run_ids=[]
    for col in header:
        if col.startswith('membrane | V (millivolt) | Run #'):
            run_ids.append(int(col.rsplit('#',1)[1]))
    if not run_ids: raise ValueError('CSV: falta membrane | V (millivolt) | Run #X')
    return t, vals, header, sorted(set(run_ids))


def column(h, a, key:str, run:int):
    q=f'{key} | Run #{run}'
    if q not in h: raise ValueError('CSV: falta '+q)
    return a[:,h.index(q)]


def run_hhlab(rhs,t,init):
    # El protocolo es una función a trozos. Integrar por intervalos evita
    # atravesar discontinuidades 10 y 10.5 ms sin control.
    intervals=[(float(t[0]),10.,0.),(10.,10.5,20.),(10.5,float(t[-1]),0.)]
    out=np.zeros((4,len(t)))
    y0=np.asarray(init,float)
    for j,(start,end,I) in enumerate(intervals):
        if end <= start: raise ValueError('El CSV debe cubrir de 0 a 45 ms y ambos cambios de estímulo')
        sol=solve_ivp(rhs,(start,end),y0,args=(I,),method='DOP853',
                      rtol=2e-11,atol=1e-13,max_step=0.05,dense_output=True)
        if not sol.success: raise RuntimeError(sol.message)
        sel=(t>=start)&(t<end if j<2 else t<=end)
        out[:,sel]=sol.sol(t[sel])
        y0=sol.y[:,-1]
    return out


def write_csv(dest,keys,rows):
    dest.parent.mkdir(parents=True,exist_ok=True)
    with dest.open('w',newline='',encoding='utf-8') as f:
        writer=csv.writer(f);writer.writerow(keys);writer.writerows(rows)


def main():
    ap=argparse.ArgumentParser(description='Validación OpenCOR contra HHlab con protocolo idéntico')
    ap.add_argument('--csv',type=Path,default=DEFAULT_CSV)
    ap.add_argument('--cellml',type=Path,default=DEFAULT_CELLML)
    ap.add_argument('--out',type=Path,default=BASE/'resultados'/'physiome_externo')
    ap.add_argument('--run',type=int,default=1,help='Ejecución para graficar (por defecto #1). Se verifican todas.')
    ap.add_argument('--hhlab-root',type=Path,default=BASE,help='Ruta de HHlab (debe contener src/hh_fourier.py)')
    args=ap.parse_args()
    args.csv=args.csv.resolve();args.cellml=args.cellml.resolve();args.out=args.out.resolve()
    hhlab_mod=(args.hhlab_root/'src'/'hh_fourier.py').resolve()
    if not hhlab_mod.is_file(): raise FileNotFoundError('Falta módulo principal HHlab: '+str(hhlab_mod))
    sys.path.insert(0,str(hhlab_mod.parent))
    from hh_fourier import hh_rhs

    model=inspect_cellml(args.cellml)
    t,a,h,runs=load_opencor(args.csv)
    if args.run not in runs: raise ValueError(f'Run #{args.run} no disponible. Hay: {runs}')
    if not np.isclose(t[0],0,atol=1e-12) or not np.isclose(t[-1],45,atol=1e-8):
        raise ValueError('Se esperaba rango [0,45] ms; revisa la exportación.')
    args.out.mkdir(parents=True,exist_ok=True)

    rows=[]; selected=None
    for run in runs:
        Vr=column(h,a,'membrane | V (millivolt)',run)
        I=column(h,a,'membrane | i_Stim (microA_per_cm2)',run)
        m=column(h,a,'sodium_channel_m_gate | m (dimensionless)',run)
        hgate=column(h,a,'sodium_channel_h_gate | h (dimensionless)',run)
        n=column(h,a,'potassium_channel_n_gate | n (dimensionless)',run)
        # CSV utilizado SOLO como referencia inicial cuando Run>1.
        # El Run #1 se valida contra el propio .cellml inicial.
        initial=[-65.-Vr[0],m[0],hgate[0],n[0]]
        if run==1 and not np.allclose(initial,[-65, .05,.6,.325],atol=1e-10):
            raise ValueError(f'Run #1 no coincide con condiciones iniciales CellML: {initial}')
        # Tiempos en los que el estímulo debe ser -20; en los extremos se tolera
        # el valor de cualquiera de los lados debido a discontinuidad.
        mask=(t>10.0)&(t<10.5)
        other=(t<10.0)|(t>10.5)
        if np.max(np.abs(I[mask]+20))>1e-9 or np.max(abs(I[other]))>1e-9:
            raise ValueError(f'Estímulo del Run #{run} no coincide con CellML')
        yp=run_hhlab(hh_rhs,t,initial)
        Vref=-65.-Vr;Vpred=yp[0];diff=Vpred-Vref
        metrics={'run':run,'muestras':int(len(t)),'rmse_mV':float(np.sqrt(np.mean(diff**2))),
                 'max_abs_mV':float(np.max(abs(diff))),
                 'bias_mV':float(np.mean(diff)),
                 'Vmax_opencor_convertido_mV':float(np.max(Vref)),
                 'Vmax_HHlab_mV':float(np.max(Vpred)),
                 't_pico_opencor_ms':float(t[np.argmax(Vref)]),
                 't_pico_HHlab_ms':float(t[np.argmax(Vpred)]),
                 'Vmin_opencor_convertido_mV':float(np.min(Vref)),
                 'Vmin_HHlab_mV':float(np.min(Vpred)),
                 'm_rmse':float(np.sqrt(np.mean((yp[1]-m)**2))),
                 'h_rmse':float(np.sqrt(np.mean((yp[2]-hgate)**2))),
                 'n_rmse':float(np.sqrt(np.mean((yp[3]-n)**2)))}
        rows.append(metrics)
        if run==args.run:
            selected=(Vref,Vpred,diff,yp,m,hgate,n,I)
    fields=list(rows[0]); write_csv(args.out/'metricas_por_ejecucion.csv',fields,[[r[k] for k in fields] for r in rows])
    Vref,Vpred,err,yp,m,hgate,n,I=selected
    write_csv(args.out/f'comparacion_run_{args.run:02d}.csv',
              ['tiempo_ms','V_opencor_original_mV','V_opencor_convertido_mV','V_HHlab_mV','error_HHlab_menos_OpenCOR_mV',
               'm_opencor','m_HHlab','h_opencor','h_HHlab','n_opencor','n_HHlab'],
              zip(t,-65-Vref,Vref,Vpred,err,m,yp[1],hgate,yp[2],n,yp[3]))

    fig,ax=plt.subplots(figsize=(10,5.4))
    ax.plot(t,Vref,lw=2,label='OpenCOR (convertido a -65 - V_original)')
    ax.plot(t,Vpred,lw=1.3,ls='--',label='HHlab Python (mismo protocolo)')
    ax.set(xlabel='Tiempo (ms)',ylabel='Voltaje transmembrana (mV)',
           title=f'Validación externa OpenCOR vs HHlab — Run #{args.run}')
    ax.grid(alpha=.3);ax.legend();fig.tight_layout()
    fig.savefig(args.out/'01_comparacion_voltaje.png',dpi=175);plt.close(fig)

    fig,ax=plt.subplots(figsize=(10,4.4));ax.plot(t,err,lw=1)
    ax.axhline(0,lw=.8,color='gray')
    ax.set(xlabel='Tiempo (ms)',ylabel='Diferencia (mV)',
           title=f'Error de voltaje: HHlab − OpenCOR — Run #{args.run}')
    ax.grid(alpha=.3);fig.tight_layout();fig.savefig(args.out/'02_error_temporal.png',dpi=175);plt.close(fig)

    fig,ax=plt.subplots(figsize=(10,5.2))
    for k,reference,prediction in zip(('m','h','n'),(m,hgate,n),(yp[1],yp[2],yp[3])):
        ax.plot(t,reference,lw=1.4,label=f'{k} OpenCOR')
        ax.plot(t,prediction,lw=1.0,ls='--',label=f'{k} HHlab')
    ax.set(xlabel='Tiempo (ms)',ylabel='Variable de compuerta (0 a 1)',
           title=f'Comparación de compuertas iónicas — Run #{args.run}')
    ax.grid(alpha=.3);ax.legend(ncol=3,fontsize=8);fig.tight_layout()
    fig.savefig(args.out/'03_compuertas.png',dpi=175);plt.close(fig)

    fig,ax=plt.subplots(figsize=(9.2,4.7))
    ax.plot([r['run'] for r in rows],[r['rmse_mV']*1000 for r in rows],marker='o')
    ax.set(xlabel='Número de ejecución exportada (Run #)',ylabel='RMSE (µV)',
           title='Concordancia de OpenCOR y HHlab por ejecución')
    ax.set_xticks([r['run'] for r in rows]);ax.grid(alpha=.3);fig.tight_layout()
    fig.savefig(args.out/'04_rmse_por_ejecucion.png',dpi=175);plt.close(fig)

    report={'estado':'comparación externa ejecutada',
            'nota':'Se ha comparado HHlab con el CSV exportado de OpenCOR, no con un segundo código Python',
            'modelo':'Modelos equivalentes para el experimento del tutorial de OpenCOR, no idéntico al estímulo constante del informe UTP.',
            'modelo_cellml':model,
            'convencion':'V_HHlab = -65 - V_CellML (mV); I_HHlab = -i_Stim_CellML',
            'muestras_por_run':len(t),'tiempo_ms_inicio':float(t[0]),'tiempo_ms_final':float(t[-1]),
            'resolucion_temporal_ms':float(np.median(np.diff(t))),
            'n_ejecuciones_csv':len(runs),'run_graficado':args.run,
            'metricas_run_elegido':next(r for r in rows if r['run']==args.run),
            'resumen_13_ejecuciones':{
                'rmse_max_mV':float(max(r['rmse_mV'] for r in rows)),
                'rmse_mediana_mV':float(np.median([r['rmse_mV'] for r in rows])),
                'max_abs_global_mV':float(max(r['max_abs_mV'] for r in rows))},
            'archivos':{'cellml':str(args.cellml),'cellml_sha256':sha256(args.cellml),
                        'csv_opencor':str(args.csv),'csv_opencor_sha256':sha256(args.csv),
                        'codigo_hhlab':str(hhlab_mod),'codigo_hhlab_sha256':sha256(hhlab_mod)},
            'advertencias':[
                'Run #2–13 pueden iniciar con el estado final previo: no son réplicas experimentales independientes.',
                'La comparación utiliza el protocolo de pulso del CellML, no la corriente constante de 10 µA/cm² del informe.',
                'Este modelo tutorial no está automáticamente identificado como el experimento concreto e/672 de Physiome.',
                'Concordancia numérica no equivale a validación fisiológica o clínica.'
            ]}
    (args.out/'reporte_comparacion.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps({'estado':'OK','run':args.run,'n_runs':len(runs),
        'rmse_mV':report['metricas_run_elegido']['rmse_mV'],
        'max_abs_mV':report['metricas_run_elegido']['max_abs_mV'],
        't_pico_ms':report['metricas_run_elegido']['t_pico_HHlab_ms'],
        'salida':str(args.out)},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
