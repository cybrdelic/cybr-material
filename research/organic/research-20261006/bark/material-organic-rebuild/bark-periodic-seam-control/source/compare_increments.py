import json
from periodic_sector import ROOT
a=json.loads((ROOT/'receipts/sector48.json').read_text());b=json.loads((ROOT/'receipts/sector48_32.json').read_text())
fa=a['trace'][-1];fb=b['trace'][-1]
result=dict(status='PASS',load=1.,intervals=[16,32],
    stored_energy_relative_change=abs(fa['stored_energy_J']-fb['stored_energy_J'])/fb['stored_energy_J'],
    imposed_work_relative_change=abs(fa['eigenstrain_work_J']-fb['eigenstrain_work_J'])/fb['eigenstrain_work_J'],
    work_ledger_relative_error=[abs(fa['energy_work_error_J'])/fa['stored_energy_J'],abs(fb['energy_work_error_J'])/fb['stored_energy_J']],
    finer_stored_energy_J=fb['stored_energy_J'],finer_work_J=fb['eigenstrain_work_J'],
    scope='Quasistatic load-increment comparison with unchanged material law; no physical rate effect.')
assert result['stored_energy_relative_change']<1e-7 and result['imposed_work_relative_change']<.001
assert result['work_ledger_relative_error'][1]<result['work_ledger_relative_error'][0]
(ROOT/'receipts/load_increment_comparison.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
