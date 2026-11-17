"""Descarga EXPLÍCITA de grabaciones reales Allen Cell Types mediante AllenSDK.

Este script depende de Internet y AllenSDK; NO genera datos falsos.
Registra ID, sweep, especie (si hay metadata), unidades, fs, SHA-256 NWB y URL.
Usa varios barridos de estímulo Long Square y clips de [-4,+8) ms del primer spike.
El CSV generado conserva muestras de la adquisición nativa (sin interpolación).
"""
import argparse
import csv
import hashlib
import json
import traceback
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def sha256_file(path):
    digest=hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(2**20),b''):digest.update(chunk)
    return digest.hexdigest()

def export(specimens,out,max_sweeps=2,window_before_ms=4.,window_after_ms=8.):
    try: from allensdk.core.cell_types_cache import CellTypesCache
    except Exception as e: raise SystemExit('AllenSDK no disponible en este Python. Usar entorno separado Python 3.10 (véase docs/ALLEN.md). '+str(e))
    out.mkdir(parents=True,exist_ok=True)
    cache=CellTypesCache(manifest_file=str(out/'manifest.json'))
    info=[]; failures=[]
    for specimen_id in specimens:
        try:
            # Consultar el espécimen explícito, no adjudicar especie por deducción.
            metas=cache.get_cells()
            met=next((x for x in metas if x.get('id')==specimen_id),{})
            ds=cache.get_ephys_data(specimen_id)
            sw=cache.get_ephys_sweeps(specimen_id)
            candidates=[x for x in sw if 'Long Square' in x.get('stimulus_name','')]
            selected=0
            nwb=out/f'specimen_{specimen_id}'/'ephys.nwb'
            checksum=sha256_file(nwb) if nwb.exists() else None
            for item in sorted(candidates,key=lambda x:x['sweep_number']):
                if selected>=max_sweeps:break
                number=int(item['sweep_number'])
                spikes=np.asarray(ds.get_spike_times(number))
                if len(spikes)==0:continue
                raw=ds.get_sweep(number)
                fs=float(raw['sampling_rate']);response=np.asarray(raw['response'])
                stimulus=np.asarray(raw['stimulus'])
                if fs<=1000 or len(response)!=len(stimulus):continue
                spike_time=float(spikes[0]);start=max(0,int(np.floor((spike_time-window_before_ms/1000)*fs)))
                stop=min(len(response),int(np.ceil((spike_time+window_after_ms/1000)*fs)))
                idx_range=raw.get('index_range', [0,len(response)-1])
                start=max(start,int(idx_range[0]));stop=min(stop,int(idx_range[1])+1)
                if stop-start<100:continue
                rr=response[start:stop]
                if not np.all(np.isfinite(rr)) or np.ptp(rr)*1000<25:continue
                stem=f'allen_{specimen_id}_sweep_{number}'
                dest=out/f'{stem}.csv'
                with dest.open('w',newline='',encoding='utf8') as f:
                    w=csv.writer(f);w.writerow(['tiempo_ms','voltaje_mV','estimulo_pA'])
                    for ix in range(start,stop):
                        w.writerow([(ix/fs-spike_time)*1000,response[ix]*1000,stimulus[ix]*1e12])
                info.append({'tipo':'EXPERIMENTAL_REAL','specimen_id':specimen_id,'sweep_number':number,
                             'source':'Allen Cell Types Database',
                             'source_docs':'https://alleninstitute.github.io/AllenSDK/cell_types.html',
                             'sampling_rate_Hz':fs,'spike_time_s':spike_time,
                             'species':met.get('species','NO_CONSIGNADO'),
                             'stimulus_name':item.get('stimulus_name'),
                             'nwb_path':str(nwb) if nwb.exists() else 'AllenSDK cache: revisar manifest.json',
                             'nwb_sha256':checksum,'csv':str(dest.relative_to(out)),
                             'num_muestras':stop-start,'conversion':'V -> mV y A -> pA',
                             'cortado_ms':[-window_before_ms,window_after_ms]})
                selected+=1
            if not selected:failures.append({'specimen_id':specimen_id,'causa':'Sin sweeps Long Square con espiga utilizable'})
        except Exception as ex:
            traceback.print_exc()  # Informa la linea exacta si aparece otro problema
            failures.append({'specimen_id':specimen_id,'causa':str(ex)})
    (out/'procedencia_allen.json').write_text(json.dumps({'registros':info,'fallos':failures},indent=2,ensure_ascii=False),encoding='utf8')
    if not info:raise SystemExit(f'Sin archivos experimentales descargados. Errores: {failures}')
    print('CSV experimentales creados:',len(info),'| fallos:',len(failures))
    return info

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--specimens',type=int,nargs='+',default=[324257146,320654829],
                   help='Ejemplos oficiales Allen; puede sustituirlos por specimen IDs confirmados')
    p.add_argument('--max-sweeps',type=int,default=2)
    p.add_argument('--out',type=Path,default=ROOT/'data'/'experimental_allen')
    a=p.parse_args();export(a.specimens,a.out,a.max_sweeps)
