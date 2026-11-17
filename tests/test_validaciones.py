import csv
import json
import tempfile
import unittest
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from validacion_physiome import hh_relative_rhs,compare
from hh_fourier import hh_rhs,gates_equilibrium,fft_coefficients,reconstruct
from analisis_experimental import single_trace,read_trace,demo

class TestIntegracionCientifica(unittest.TestCase):
    def test_derivadas_dos_convenciones(self):
        for voltage in [-80,-65,-55,-40,10,30]:
            y=gates_equilibrium();y[0]=voltage
            a=np.asarray(hh_rhs(0,y,12))
            b=y.copy();b[0]+=65
            c=np.asarray(hh_relative_rhs(0,b,12))
            self.assertLess(np.max(np.abs(a-c)),1e-12)

    def test_equivalencia_numerica(self):
        with tempfile.TemporaryDirectory() as d:
            r=compare(Path(d),duration=50.0)
            self.assertFalse(r['physiome_cellml_ejecutado'])
            self.assertLess(r['max_abs_diferencia_mV'],.0002)
            self.assertTrue((Path(d)/'comparacion_convenciones.csv').is_file())

    def test_registro_demostracion_no_atribuido_a_allen(self):
        with tempfile.TemporaryDirectory() as d:
            inp=demo(Path(d))
            path=inp/'DEMO_SINTETICA_1.csv'
            self.assertTrue(path.is_file())
            t,v=read_trace(path)
            self.assertEqual(len(t),len(v))
            rows,ex=single_trace(path,'SINTETICO_DEMOSTRACION','demo1')
            self.assertTrue(all(x['tipo']=='SINTETICO_DEMOSTRACION' for x in rows))
            self.assertEqual(len(rows),200) # N=1..100, dos métodos
            self.assertEqual(len(ex['FFT_cruda_N25']),len(v))
            self.assertEqual(len(ex['FFT_residuo_con_puente_N25']),len(v))

    def test_clip_par_corrige_longitud(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'odd.csv'
            with path.open('w',newline='') as f:
                w=csv.writer(f);w.writerow(['tiempo_ms','voltaje_mV'])
                w.writerows((i*.02, -65+80*np.exp(-((i*.02-1)/.3)**2)) for i in range(201))
            t,v=read_trace(path)
            self.assertEqual(len(v),200)
            self.assertLess(np.max(np.abs(reconstruct(fft_coefficients(v),99)-v)),3)

    def test_puente_recupera_tendencia_lineal(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'trend.csv'
            n=200
            t=np.arange(n)*.02
            v=-65+2*t+4*np.sin(2*np.pi*5*np.arange(n)/n)
            with path.open('w',newline='') as f:
                w=csv.writer(f);w.writerow(['tiempo_ms','voltaje_mV']);w.writerows(zip(t,v))
            rows,_=single_trace(path,'SINTETICO_DEMOSTRACION','lin')
            rmse_raw=next(x['RMSE_mV'] for x in rows if x['metodo']=='FFT_cruda' and x['N']==5)
            rmse_adj=next(x['RMSE_mV'] for x in rows if x['metodo']=='FFT_residuo_con_puente' and x['N']==5)
            self.assertLess(rmse_adj,rmse_raw)

if __name__=='__main__': unittest.main()
