"""Extrae trazas *reales* de Allen Cell Types Database mediante AllenSDK.

Requiere conexión a Internet y ALLENSDK en un entorno compatible. La ejecución
local del proyecto NO necesita AllenSDK. Fuente de los datos: Allen Institute,
Cell Types Database. Los datos de ratón/humano no validan numéricamente un modelo
de axón de calamar sin ajuste de parámetros y protocolo.

Uso opcional:
    pip install allensdk
    python src/obtener_datos_allen.py --specimen-id 464212183
"""
from pathlib import Path
import argparse
import csv
import json
import numpy as np


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--specimen-id',type=int,default=464212183)
    parser.add_argument('--out',type=Path,default=Path(__file__).resolve().parents[1]/'data'/'experimental_allen')
    args=parser.parse_args()
    try:
        from allensdk.core.cell_types_cache import CellTypesCache
    except ImportError:
        raise SystemExit('Falta allensdk: instalarlo en un venv compatible. Ver docs/BASES_DE_DATOS.md')
    args.out.mkdir(parents=True,exist_ok=True)
    cache=CellTypesCache(manifest_file=str(args.out/'manifest.json'))
    ds=cache.get_ephys_data(args.specimen_id)
    sweeps=cache.get_ephys_sweeps(args.specimen_id)
    candidates=[s for s in sweeps if 'Long Square' in s.get('stimulus_name','')]
    if not candidates:
        raise RuntimeError('No hay sweeps Long Square. Inspeccione metadatos Allen para este specimen.')
    selection=None
    for item in candidates:
        number=item['sweep_number']
        spikes=ds.get_spike_times(number)
        if len(spikes)>0:
            selection=(number,float(spikes[0]));break
    if selection is None:
        raise RuntimeError('Ningún sweep Long Square seleccionado presenta espigas.')
    number, spike_time=selection
    raw=ds.get_sweep(number)
    fs=float(raw['sampling_rate'])
    response=np.asarray(raw['response'])  # unidades Voltios: convertir a mV
    stimulus=np.asarray(raw['stimulus']) # unidades Amperios: convertir a pA
    ti=np.asarray(raw.get('index_range',[0,len(response)-1]))
    left=max(int(ti[0]),int((spike_time-.004)*fs))
    right=min(int(ti[1])+1,int((spike_time+.011)*fs))
    if right-left<100:raise RuntimeError('Ventana demasiado pequeña para el barrido seleccionado.')
    out=args.out/f'allen_{args.specimen_id}_sweep_{number}.csv'
    with out.open('w',newline='',encoding='utf8') as f:
        w=csv.writer(f)
        w.writerow(['tiempo_ms','voltaje_mV','estimulo_pA'])
        w.writerows(zip(np.arange(left,right)*1000/fs, response[left:right]*1000,
                        stimulus[left:right]*1e12))
    info={'origen':'Allen Cell Types Database, patch-clamp intracelular',
          'specimen_id':args.specimen_id,'sweep_number':number,
          'sampling_rate_Hz':fs,'spike_time_s':spike_time,
          'tipo_datos':'experimental real (NO simulado)',
          'especie':'consultar metadatos Allen del espécimen; no asumir humano',
          'url_fuente':'https://celltypes.brain-map.org/',
          'no_equivalencia':'No calibrado al calamar: no comparar voltajes punto a punto como validación HH 1952'}
    (args.out/'procedencia.json').write_text(json.dumps(info,indent=2,ensure_ascii=False),encoding='utf8')
    print('CSV experimental descargado y generado:',out)

if __name__=='__main__':main()
