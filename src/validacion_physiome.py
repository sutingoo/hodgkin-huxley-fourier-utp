"""Validación independiente EN CONVENCIONES: HH moderno vs relativo al reposo.

IMPORTANTE: NO ejecuta el CellML de Physiome. Verifica la equivalencia
algebraica y numérica de dos implementaciones con unidades explícitas.
Con --opencor-csv realiza además comparación externa *si* se exporta una
simulación OpenCOR con protocolo y parámetros coincidentes.
"""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import solve_ivp
from hh_fourier import gates_equilibrium, simulate

ROOT=Path(__file__).resolve().parents[1]

def _trap_classic(a:float, scale:float)->float:
    # a/(exp(a/scale)-1) ; límite scale. Forma distinta al modelo moderno.
    z=a/scale
    if abs(z)<1e-6:
        return scale*(1-z/2+z*z/12)
    return float(a/np.expm1(z))

def hh_relative_rhs(t, y, I):
    """v=Vmodern+65; E_Na=115,E_K=-12,E_L=10.613 (mV relativos)."""
    v,m,h,n=y
    am=0.1*_trap_classic(25-v,10)
    bm=4*np.exp(-v/18)
    ah=0.07*np.exp(-v/20)
    bh=1/(1+np.exp((30-v)/10))
    an=0.01*_trap_classic(10-v,10)
    bn=0.125*np.exp(-v/80)
    ina=120*m**3*h*(v-115)
    ik=36*n**4*(v+12)
    il=.3*(v-10.613)
    return [(I-ina-ik-il), am*(1-m)-bm*m, ah*(1-h)-bh*h, an*(1-n)-bn*n]

def compare(output:Path, I=10.0, duration=100.0, opencor_csv=None, time_col='time_ms',
            voltage_col='voltage_mV', reference_offset=-65.0):
    output.mkdir(parents=True,exist_ok=True)
    modern=simulate(current=I,duration_ms=duration,method='DOP853',rtol=2e-10,atol=1e-12)
    y0=gates_equilibrium(); y0[0]+=65
    relative=solve_ivp(hh_relative_rhs,(0,duration),y0,args=(I,),
                       method='RK45',rtol=2e-10,atol=1e-12,max_step=.12,dense_output=True)
    if not relative.success: raise RuntimeError(relative.message)
    t=np.arange(0,duration,.05)
    a=modern.sol(t);b=relative.sol(t);b[0]-=65
    err=b[0]-a[0]
    with (output/'comparacion_convenciones.csv').open('w',newline='',encoding='utf8') as f:
        w=csv.writer(f);w.writerow(['tiempo_ms','V_modelo_moderno_mV','V_modelo_relativo_convertido_mV','error_mV'])
        w.writerows(zip(t,a[0],b[0],err))
    report={'tipo':'validacion interna: implementaciones independientes en convención de potencial',
            'physiome_cellml_ejecutado':False,
            'duracion_ms':duration,'I_ext_uA_cm2':I,
            'max_abs_diferencia_mV':float(np.max(np.abs(err))),
            'rmse_diferencia_mV':float(np.sqrt(np.mean(err**2))),
            'max_diferencia_compuertas':float(np.max(np.abs(a[1:]-b[1:]))),
            'advertencia':'No es una ejecución de Physiome/OpenCOR; hay que exportar el modelo oficial y armonizar el estímulo.'}
    if opencor_csv:
        dat=np.genfromtxt(opencor_csv,delimiter=',',names=True,encoding='utf8')
        if time_col not in dat.dtype.names or voltage_col not in dat.dtype.names:
            raise ValueError(f'Columnas del CSV: {dat.dtype.names}. Se esperaban {time_col}, {voltage_col}')
        tt=np.asarray(dat[time_col],dtype=float)
        vv=np.asarray(dat[voltage_col],dtype=float)+reference_offset
        valid=np.isfinite(tt)&np.isfinite(vv)&(tt>=0)&(tt<=duration)
        tt=tt[valid];vv=vv[valid]
        if len(tt)<50 or np.any(np.diff(tt)<=0): raise ValueError('Tiempos exportados incorrectos')
        d=vv-modern.sol(tt)[0]
        report['opencor_comparacion']={'csv':str(opencor_csv),'columna':voltage_col,
            'offset_convertir_a_absoluto_mV':reference_offset,
            'muestras_comparadas':len(tt),'rmse_mV':float(np.sqrt(np.mean(d*d))),
            'max_abs_mV':float(np.max(np.abs(d)))}
        report['physiome_cellml_ejecutado']='CSV externo aportado; debe comprobarse su procedencia y protocolo'
    (output/'reporte_physiome.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf8')
    return report

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--out',type=Path,default=ROOT/'resultados'/'validacion')
    p.add_argument('--opencor-csv',type=Path)
    p.add_argument('--time-col',default='time_ms')
    p.add_argument('--voltage-col',default='voltage_mV')
    p.add_argument('--reference-offset',type=float,default=-65.0,
                   help='Si Physiome exporta V relativo a reposo, -65; si es absoluto, 0.')
    p.add_argument('--current',type=float,default=10.)
    args=p.parse_args()
    print(json.dumps(compare(args.out,args.current,opencor_csv=args.opencor_csv,
          time_col=args.time_col,voltage_col=args.voltage_col,
          reference_offset=args.reference_offset),indent=2,ensure_ascii=False))
