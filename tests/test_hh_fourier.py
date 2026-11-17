import sys
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from hh_fourier import (vtrap, rates, gates_equilibrium, hh_rhs, simulate, isolated_cycle,
                        fft_coefficients, reconstruct, ac_energy_fraction, metrics)

class TestHodgkinHuxleyFourier(unittest.TestCase):
    def test_singularidades_removibles(self):
        self.assertAlmostEqual(vtrap(0.0,10.0),10.0,places=10)
        self.assertAlmostEqual(rates(-40)[0],1.0,places=10)
        self.assertAlmostEqual(rates(-55)[4],.1,places=10)

    def test_equilibrio_compuertas(self):
        s = gates_equilibrium()
        self.assertTrue(np.all((s[1:]>0)&(s[1:]<1)))
        derivative = hh_rhs(0,s,0)
        self.assertLess(np.max(np.abs(derivative[1:])),1e-12)
        self.assertLess(abs(derivative[0]),0.02)

    def test_fourier_senoides_exacto(self):
        t=np.arange(4096)/4096
        v=-60+15*np.cos(2*np.pi*3*t)+3*np.sin(2*np.pi*7*t)
        c=fft_coefficients(v)
        self.assertLess(np.max(np.abs(reconstruct(c,7)-v)),1e-10)
        self.assertAlmostEqual(np.sqrt(np.mean((reconstruct(c,3)-v)**2)),3/np.sqrt(2),places=10)
        self.assertAlmostEqual(ac_energy_fraction(c,7),1.0,places=12)

    def test_parseval_y_error_residual(self):
        t=np.arange(4096)/4096
        v=-60+np.cos(2*np.pi*t)+2*np.sin(2*np.pi*12*t)
        c=fft_coefficients(v)
        # Parseval de la componente AC, consistente con el RMSE del error.
        for N in [1,5,10,12]:
            r=reconstruct(c,N)
            e=1-ac_energy_fraction(c,N)
            acpow=np.mean((v-np.mean(v))**2)
            self.assertAlmostEqual(np.mean((r-v)**2),e*acpow,places=9)

    def test_disparo_y_fidelidad(self):
        sol=simulate(current=10,duration_ms=160)
        t,y,T,p=isolated_cycle(sol)
        c=fft_coefficients(y[0]); r25=reconstruct(c,25)
        met=metrics(t,y[0],r25,c,25,T)
        self.assertTrue(14<T<15)
        self.assertTrue(25<np.max(y[0])<40)
        self.assertTrue(-80<np.min(y[0])<-70)
        self.assertLess(met['RMSE_mV'],0.8)
        self.assertGreater(met['energia_AC_retenida_pct'],99.9)

if __name__ == '__main__':unittest.main()
