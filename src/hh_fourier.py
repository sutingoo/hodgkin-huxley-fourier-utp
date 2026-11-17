"""Hodgkin-Huxley clásico (1952) + análisis Fourier, reproducible.

Unidad tiempo ms; V mV; corrientes uA/cm2; conductancias mS/cm2.
Se usa la parametrización trasladada a Vreposo=-65 mV.
No modela una neurona humana ni un electrodo extracelular.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
from scipy.optimize import minimize_scalar

CM = 1.0
G_NA, G_K, G_L = 120.0, 36.0, 0.3
E_NA, E_K, E_L = 50.0, -77.0, -54.387
V_REST = -65.0
N_SAMPLES = 4096
N_MAX = 100
CANDIDATES = [1, 3, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100]


def vtrap(x: float, scale: float) -> float:
    """Evalúa x/(1-exp(-x/scale)) sin singularidad 0/0 en x=0."""
    z = x / scale
    if abs(z) < 1e-6:
        return scale * (1 + z / 2 + z*z / 12)
    return float(x / (-np.expm1(-z)))


def rates(V: float) -> tuple[float, float, float, float, float, float]:
    return (
        0.1 * vtrap(V + 40.0, 10.0),
        4.0 * np.exp(-(V + 65.0) / 18.0),
        0.07 * np.exp(-(V + 65.0) / 20.0),
        1.0 / (1.0 + np.exp(-(V + 35.0) / 10.0)),
        0.01 * vtrap(V + 55.0, 10.0),
        0.125 * np.exp(-(V + 65.0) / 80.0),
    )


def gates_equilibrium(V: float = V_REST) -> np.ndarray:
    am, bm, ah, bh, an, bn = rates(V)
    return np.array([V, am/(am+bm), ah/(ah+bh), an/(an+bn)])


def ionic_currents(y: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    V, m, h, n = y
    return G_NA*(m**3)*h*(V-E_NA), G_K*(n**4)*(V-E_K), G_L*(V-E_L)


def hh_rhs(_t: float, y: np.ndarray, current: float) -> list[float]:
    V, m, h, n = y
    am, bm, ah, bh, an, bn = rates(float(V))
    ina, ik, il = ionic_currents(y)
    return [(current-ina-ik-il)/CM,
            am*(1-m)-bm*m, ah*(1-h)-bh*h, an*(1-n)-bn*n]


def simulate(current: float = 10.0, duration_ms: float = 250.0,
             method: str = 'DOP853', rtol: float = 2e-10, atol: float = 1e-12):
    sol = solve_ivp(hh_rhs, (0, duration_ms), gates_equilibrium(),
                    args=(current,), method=method, rtol=rtol, atol=atol,
                    max_step=0.12, dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def sample_trace(sol, current: float, sample_dt: float = 0.05):
    t = np.arange(0, sol.t[-1], sample_dt)
    y = sol.sol(t)
    ina, ik, il = ionic_currents(y)
    return {'tiempo_ms': t, 'voltaje_mV': y[0], 'm': y[1], 'h': y[2], 'n': y[3],
            'I_ext_uA_cm2': np.full(len(t), current),
            'I_Na_uA_cm2': ina, 'I_K_uA_cm2': ik, 'I_L_uA_cm2': il}


def isolated_cycle(sol, n_samples: int = N_SAMPLES):
    # Más de 2 ciclos completos en el tramo estacionario, selección automatizada.
    t_probe = np.arange(max(30.0, sol.t[-1]-125), sol.t[-1], 0.025)
    v_probe = sol.sol(t_probe)[0]
    inds, _ = find_peaks(v_probe, height=0, prominence=35, distance=int(8/0.025))
    if len(inds) < 3:
        raise ValueError('No se detectaron al menos tres espigas periódicas; ajuste I_ext.')
    rough_peaks = t_probe[inds]
    peaks = np.array([minimize_scalar(lambda z: -float(sol.sol(z)[0]),
                                     bounds=(float(x-0.04), float(x+0.04)),
                                     method="bounded", options={"xatol":1e-11}).x
                      for x in rough_peaks])
    T = float(peaks[-1] - peaks[-2])
    if not (7.0 < T < 30.0):
        raise ValueError(f'Periodo inusual: {T} ms')
    # Colocar el pico alrededor de t_local=3 ms.
    t0 = float(peaks[-2] - 3.0)
    t = np.arange(n_samples) * T / n_samples
    y = sol.sol(t0 + t)
    if np.max(np.abs(sol.sol(np.array([t0, t0+T]))[:, 0] -
                         sol.sol(np.array([t0, t0+T]))[:, 1])) > 0.25:
        raise ValueError('La onda aún no es suficientemente periódica para DFT de un ciclo.')
    return t, y, T, peaks


def fft_coefficients(v: np.ndarray) -> np.ndarray:
    return np.fft.rfft(v) / len(v)


def reconstruct(c: np.ndarray, N: int, lanczos: bool = False) -> np.ndarray:
    n = (len(c)-1)*2
    out = np.zeros_like(c)
    out[:N+1] = c[:N+1]
    if lanczos:
        k = np.arange(1, N+1)
        out[1:N+1] *= np.sinc(k/(N+1))
    return np.fft.irfft(out * n, n=n)


def ac_energy_fraction(c: np.ndarray, N: int) -> float:
    # Coeficientes de 0 y Nyquist no pertenecen a la energía AC acumulada por armónicos.
    weights = np.full(len(c), 2.0)
    weights[0] = 0.0
    weights[-1] = 1.0
    e = weights * (np.abs(c)**2)
    return float(e[1:N+1].sum() / e[1:].sum())


def crossing_time(t: np.ndarray, v: np.ndarray, threshold: float=-20.0):
    eligible = np.flatnonzero((v[:-1] < threshold) & (v[1:] >= threshold))
    if len(eligible)==0:
        return float('nan')
    i = int(eligible[0])
    return float(t[i]+(threshold-v[i])*(t[i+1]-t[i])/(v[i+1]-v[i]))


def metrics(t: np.ndarray, v: np.ndarray, rec: np.ndarray, coeffs, N: int, T:float):
    err = rec-v
    vpp = float(np.ptp(v))
    rmse = float(np.sqrt(np.mean(err**2)))
    peak_err = abs(float(np.max(rec)-np.max(v)))
    max_err = float(np.max(np.abs(err)))
    overshoot = max(0.0, float(np.max(rec)-np.max(v)), float(np.min(v)-np.min(rec)))
    return {'N': N, 'f_max_Hz': 1000.0*N/T,
            'RMSE_mV': rmse, 'NRMSE_pct_Vpp': rmse/vpp*100,
            'error_pico_mV': peak_err,
            'error_max_mV': max_err, 'error_max_pct_Vpp': max_err/vpp*100,
            'sobreimpulso_pct_Vpp': overshoot/vpp*100,
            'energia_AC_retenida_pct': ac_energy_fraction(coeffs,N)*100,
            'error_cruce_umbral_ms': abs(crossing_time(t,rec)-crossing_time(t,v))}


def write_csv(path: Path, data: dict):
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(data)
    rows = zip(*(np.atleast_1d(data[k]) for k in keys))
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f); writer.writerow(keys)
        writer.writerows(rows)


def write_rows(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)


def charts(out: Path, trace, t, y, T, c, rows, recon):
    figs = out/'resultados'/'figuras'; figs.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size':10, 'axes.spines.top':False,
                         'axes.spines.right':False, 'figure.dpi':140})
    fig, ax=plt.subplots(figsize=(11,4))
    mask = trace['tiempo_ms']<=80
    ax.plot(trace['tiempo_ms'][mask], trace['voltaje_mV'][mask], lw=1.25)
    ax.set(xlabel='Tiempo (ms)', ylabel='Potencial de membrana (mV)',
           title='Modelo Hodgkin–Huxley | corriente constante 10 µA/cm²')
    ax.grid(alpha=.22); fig.tight_layout(); fig.savefig(figs/'01_tren_potenciales.png'); plt.close(fig)

    fig, ax=plt.subplots(figsize=(11,4.5))
    ax.plot(t, y[0], lw=2, label='V(t)')
    ax.plot(t, -65+95*y[1], lw=1.25, label='m(t) (escala visual)')
    ax.plot(t, -65+95*y[2], lw=1.25, label='h(t) (escala visual)')
    ax.plot(t, -65+95*y[3], lw=1.25, label='n(t) (escala visual)')
    ax.set(xlabel='Tiempo del ciclo (ms)', ylabel='Voltaje (mV) / compuertas reescaladas',
           title=f'Potencial y compuertas — T = {T:.4f} ms')
    ax.legend(ncols=4, fontsize=8); ax.grid(alpha=.2); fig.tight_layout()
    fig.savefig(figs/'02_potencial_y_compuertas.png'); plt.close(fig)

    fig,ax=plt.subplots(figsize=(11,4.5))
    ax.plot(t, y[0], color='black', lw=2, label='Referencia HH')
    for N in [3,5,10,25,50]:
        ax.plot(t, recon[N], label=f'N={N}', alpha=.8, lw=1)
    ax.set(xlabel='Tiempo del ciclo (ms)', ylabel='Potencial de membrana (mV)',
           title='Reconstrucción mediante suma parcial de Fourier')
    ax.legend(ncols=3); ax.grid(alpha=.2); fig.tight_layout()
    fig.savefig(figs/'03_reconstrucciones.png'); plt.close(fig)

    fig,ax=plt.subplots(figsize=(11,4.5))
    mask=(t>1.5)&(t<5.5)
    ax.plot(t[mask], y[0][mask], color='black', lw=2, label='Referencia HH')
    for N in [5,10,25,50]:
        ax.plot(t[mask], recon[N][mask], label=f'N={N}', lw=1.25)
    ax.set(xlabel='Tiempo (ms)', ylabel='Voltaje (mV)',title='Detalle de espiga y error por truncamiento')
    ax.legend(); ax.grid(alpha=.2); fig.tight_layout(); fig.savefig(figs/'04_detalle_espiga.png'); plt.close(fig)

    ff=np.arange(len(c))*1000/T
    fig,ax=plt.subplots(figsize=(11,4))
    ax.stem(ff[1:101], 2*np.abs(c[1:101]), basefmt=' ', markerfmt='.', linefmt='-')
    ax.set(xlabel='Frecuencia (Hz)',ylabel='Amplitud armónica (mV)',
           title='Espectro de magnitud (primeros 100 armónicos)'); ax.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(figs/'05_espectro.png'); plt.close(fig)

    fig,ax=plt.subplots(figsize=(11,4.5))
    ns=np.array([r['N'] for r in rows]); errors=np.array([r['RMSE_mV'] for r in rows]);
    ax.semilogy(ns, errors, 'o-', label='RMSE (mV)')
    ax.axhline(1.0, color='gray', ls='--', label='Criterio RMSE 1 mV')
    ax.set(xlabel='Número de armónicos N', ylabel='RMSE (mV)', title='Convergencia del error')
    ax.grid(alpha=.25); ax.legend(); fig.tight_layout(); fig.savefig(figs/'06_convergencia.png'); plt.close(fig)

    fig,ax=plt.subplots(figsize=(11,4.5))
    ns_all=np.arange(1,101); e=np.array([ac_energy_fraction(c,int(N))*100 for N in ns_all]);
    ax.plot(ns_all,e, lw=2); ax.axhline(99.9,color='gray',ls='--',label='Criterio 99.9 %')
    ax.set(xlabel='Número de armónicos N',ylabel='Energía AC retenida (%)',
           title='Teorema de Parseval: energía alterna retenida')
    ax.grid(alpha=.25); ax.legend(); fig.tight_layout(); fig.savefig(figs/'07_energia.png'); plt.close(fig)

    fig,ax=plt.subplots(figsize=(11,4.5))
    ax.plot(t,y[0],lw=2,label='Referencia HH (limpia)')
    ax.plot(t,recon[10],lw=1,label='Truncada N=10')
    ax.plot(t,reconstruct(c,10,lanczos=True),lw=1,label='Lanczos N=10')
    ax.set(xlabel='Tiempo (ms)',ylabel='Voltaje (mV)',title='Lanczos: reducción del rizado frente al sesgo')
    ax.legend(); ax.grid(alpha=.25); fig.tight_layout(); fig.savefig(figs/'08_lanczos.png'); plt.close(fig)


def run(output: Path, current:float=10.0, duration:float=250.0):
    output.mkdir(parents=True,exist_ok=True)
    sol=simulate(current=current,duration_ms=duration)
    trace=sample_trace(sol,current)
    write_csv(output/'data'/'simulados'/'hh_tren_completo.csv',trace)
    t,y,T,peaks=isolated_cycle(sol)
    c=fft_coefficients(y[0]); v=y[0]
    recon={N:reconstruct(c,N) for N in CANDIDATES}
    rows=[metrics(t,v,recon[N],c,N,T) for N in CANDIDATES]
    write_rows(output/'resultados'/'metricas_armonicos.csv',rows)
    write_csv(output/'data'/'simulados'/'hh_periodo.csv',{
        'tiempo_ms':t, 'voltaje_mV':v, 'm':y[1], 'h':y[2], 'n':y[3],
        'I_ext_uA_cm2':np.full(len(t),current)})
    write_csv(output/'resultados'/'reconstrucciones.csv',{
        'tiempo_ms':t, 'original_mV':v, **{f'Fourier_N{N}_mV':recon[N] for N in CANDIDATES},
        'Lanczos_N10_mV':reconstruct(c,10,lanczos=True)})
    ks=np.arange(1,N_MAX+1)
    write_csv(output/'resultados'/'espectro.csv',{
        'armonico_k':ks,'frecuencia_Hz':ks*1000/T,
        'amplitud_mV':2*np.abs(c[1:N_MAX+1]),
        'energia_AC_acumulada_pct':np.array([ac_energy_fraction(c,int(i))*100 for i in ks])})
    # Noise experiment on the same exact synthetic waveform.
    rng=np.random.default_rng(20261116)
    noise=rng.normal(0,1.0,len(v))
    cv=fft_coefficients(v+noise)
    write_csv(output/'resultados'/'ruido_sintetico.csv',{
        'tiempo_ms':t,'referencia_limpia_mV':v,
        'ruido_sintetico_mV':noise,'senal_ruidosa_mV':v+noise,
        'reconstruccion_N25_mV':reconstruct(cv,25)})
    e_min=None
    for n in range(1,N_MAX+1):
        met=metrics(t,v,reconstruct(c,n),c,n,T)
        if met['RMSE_mV'] < 1.0 and met['energia_AC_retenida_pct']>=99.9 and met['error_pico_mV']<1.0:
            e_min=n;break
    highres=simulate(current=current,duration_ms=duration,method='RK45',rtol=1e-10,atol=1e-12)
    check_times=np.linspace(duration-50,duration-2,1000)
    solver_difference=float(np.max(np.abs(sol.sol(check_times)[0]-highres.sol(check_times)[0])))
    summary={
       'estado':'reproducido por codigo, datos generados por simulacion; NO datos clinicos',
       'modelo':'HH1952 parametros clasicos, parametros referidos a Vrest=-65mV',
       'I_ext_uA_cm2':current,'duracion_simulacion_ms':duration,
       'periodo_ms':T,'frecuencia_fundamental_Hz':1000/T,
       'Vmax_mV':float(np.max(v)),'Vmin_mV':float(np.min(v)),
       'Vpp_mV':float(np.ptp(v)),'pendiente_max_mV_ms':float(np.max(np.gradient(v,t))),
       'muestras_por_ciclo':len(t),'paso_por_ciclo_ms':float(t[1]-t[0]),
       'condicion_optimo':'menor N: RMSE <1mV, energia_AC>=99.9%, error_pico<1mV',
       'N_min_criterio':e_min,
       'diferencia_max_solver_RK45_vs_DOP853_mV':solver_difference,
       'error_periodicidad_mV':float(np.max(np.abs(sol.sol(np.array([peaks[-2]-3.0,peaks[-2]-3.0+T]))[:,0]-sol.sol(np.array([peaks[-2]-3.0,peaks[-2]-3.0+T]))[:,1]))),
       'metricas_N25':next(x for x in rows if x['N']==25),
       'advertencia':'No prueba exactitud clinica, ahorro energetico de hardware ni optimizacion universal.'}
    (output/'resultados'/'resumen_validacion.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
    charts(output,trace,t,y,T,c,rows,recon)
    print(json.dumps(summary,indent=2,ensure_ascii=False))
    return summary

if __name__=='__main__':
    ap=argparse.ArgumentParser(description='Hodgkin-Huxley y Fourier (evidencia reproducible)')
    ap.add_argument('--out',default=str(Path(__file__).resolve().parents[1]))
    ap.add_argument('--current',type=float,default=10.0)
    ap.add_argument('--duration',type=float,default=250.0)
    args=ap.parse_args()
    run(Path(args.out),args.current,args.duration)
