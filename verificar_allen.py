#!/usr/bin/env python3
"""Audita la trazabilidad local de datos Allen exportados por HHlab.

No descarga datos. Comprueba integridad SHA-256 del NWB, CSV, muestreo y
concordancia con procedencia_allen.json. No certifica identidad biológica:
esto requiere contrastar metadatos con la API oficial del Allen Institute.

Uso desde ~/Downloads/HHlab:
python verificar_allen.py
python verificar_allen.py --manifest data/experimental_allen/procedencia_allen.json
"""
import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def inspect_csv(path: Path, expected_samples: int, expected_fs: float):
    outcome = {'archivo': str(path), 'existe': path.is_file(), 'comprobaciones': [], 'avisos': []}
    if not path.is_file():
        outcome['comprobaciones'].append('ERROR: CSV no encontrado')
        return outcome, False
    required = {'tiempo_ms', 'voltaje_mV', 'estimulo_pA'}
    count, t0, t_last, prev_t = 0, None, None, None
    bad_values, bad_steps = 0, 0
    dt_ref = 1000.0 / expected_fs
    with path.open(newline='', encoding='utf-8-sig') as stream:
        reader = csv.DictReader(stream)
        cols = set(reader.fieldnames or [])
        if not required <= cols:
            outcome['comprobaciones'].append('ERROR: columnas requeridas ausentes: ' + str(sorted(required-cols)))
            return outcome, False
        for row in reader:
            count += 1
            try:
                t = float(row['tiempo_ms']); v = float(row['voltaje_mV']); i = float(row['estimulo_pA'])
                if not all(map(math.isfinite, (t, v, i))):
                    bad_values += 1
                if t0 is None: t0 = t
                if prev_t is not None:
                    dt = t - prev_t
                    if not math.isclose(dt, dt_ref, rel_tol=2e-3, abs_tol=1e-7):
                        bad_steps += 1
                prev_t, t_last = t, t
            except (TypeError, ValueError, KeyError):
                bad_values += 1
    outcome['muestras_leidas'] = count
    outcome['muestras_esperadas'] = expected_samples
    outcome['t_inicial_ms'] = t0
    outcome['t_final_ms'] = t_last
    outcome['paso_esperado_ms'] = dt_ref
    ok = (count == expected_samples and not bad_values and not bad_steps)
    if count != expected_samples:
        outcome['comprobaciones'].append(f'ERROR: {count} filas frente a {expected_samples} esperadas')
    if bad_values:
        outcome['comprobaciones'].append(f'ERROR: {bad_values} filas con datos no finitos o inválidos')
    if bad_steps:
        outcome['comprobaciones'].append(f'ERROR: {bad_steps} pasos temporales incompatibles con fs')
    if ok:
        outcome['comprobaciones'].append('OK: longitud, columnas, valores y frecuencia de muestreo')
    return outcome, ok


def locate_nwb(record, base: Path):
    specified = Path(record.get('nwb_path', ''))
    alternatives = [specified, base / f"specimen_{record['specimen_id']}" / 'ephys.nwb']
    for path in alternatives:
        if path.is_file():
            return path
    return alternatives[-1]


def audit(manifest_path: Path):
    if not manifest_path.is_file():
        raise SystemExit(f'Manifiesto no encontrado: {manifest_path}')
    base = manifest_path.parent
    obj = json.loads(manifest_path.read_text(encoding='utf-8'))
    records = obj.get('registros', [])
    if not isinstance(records, list) or not records:
        raise SystemExit('Manifiesto sin lista de registros.')
    report = {'manifiesto': str(manifest_path.resolve()),
              'n_registros': len(records),
              'n_especimenes': len(set(x['specimen_id'] for x in records)),
              'especies_declaradas': dict(Counter(x.get('species', 'SIN_DATO') for x in records)),
              'nwb': [], 'csv': [], 'errores': [], 'avisos': [],
              'alcance': 'Auditoría local de integridad; no valida directamente metadatos contra Allen API.'}
    cache = {}
    ids_seen = set()
    for r in records:
        sid, sweep = r['specimen_id'], r['sweep_number']
        rid = (sid, sweep)
        if rid in ids_seen:
            report['errores'].append(f'Registro repetido: specimen {sid}, sweep {sweep}')
        ids_seen.add(rid)
        csv_path = base / Path(r['csv']).name
        co, ok = inspect_csv(csv_path, int(r['num_muestras']), float(r['sampling_rate_Hz']))
        co['specimen_id'], co['sweep_number'] = sid, sweep
        report['csv'].append(co)
        if not ok: report['errores'].append(f'CSV falló: {csv_path.name}')
        if r.get('stimulus_name') != 'Long Square':
            report['avisos'].append(f'Barrera {sid}/{sweep}: estímulo no Long Square (verificar)')
        nwb_path = locate_nwb(r, base)
        expected = r.get('nwb_sha256')
        file_key = str(nwb_path)
        if file_key not in cache:
            if not nwb_path.is_file():
                status = {'specimen_id': sid, 'ruta': str(nwb_path), 'estado': 'NO_ENCONTRADO'}
                report['errores'].append(f'NWB ausente: {nwb_path}')
            else:
                actual = sha256_file(nwb_path)
                status = {'specimen_id': sid, 'ruta': str(nwb_path), 'sha256_calculado': actual,
                          'sha256_manifest': expected,
                          'estado': 'OK' if actual == expected else 'HASH_DISTINTO'}
                if actual != expected: report['errores'].append(f'SHA-256 diferente para {nwb_path}')
            cache[file_key] = status
            report['nwb'].append(status)
        else:
            if expected != cache[file_key].get('sha256_manifest'):
                report['errores'].append(f'Diferentes SHA-256 declarados para NWB {nwb_path}')
    if obj.get('fallos'):
        report['avisos'].append(f"Manifiesto tiene {len(obj['fallos'])} fallos de descarga")
    report['estado'] = 'OK' if not report['errores'] else 'CON_OBSERVACIONES'
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path,
                   default=Path('data/experimental_allen/procedencia_allen.json'))
    p.add_argument('--out', type=Path,
                   default=Path('resultados/experimental/auditoria_procedencia.json'))
    args = p.parse_args()
    report = audit(args.manifest)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding='utf-8')
    print('Registros:', report['n_registros'], '| especímenes:', report['n_especimenes'])
    print('Archivos NWB examinados:', len(report['nwb']))
    for item in report['nwb']:
        print('  ',item['specimen_id'], item['estado'])
    print('CSV verificados:', sum(not any('ERROR' in s for s in x['comprobaciones']) for x in report['csv']), '/', len(report['csv']))
    print('Estado:', report['estado'])
    for s in report['errores']: print('ERROR:', s)
    print('Reporte guardado:', args.out)
    raise SystemExit(0 if report['estado'] == 'OK' else 2)

if __name__ == '__main__':
    main()
