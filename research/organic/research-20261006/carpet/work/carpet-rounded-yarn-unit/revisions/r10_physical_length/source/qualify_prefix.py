"""Read-only native gate of rejected actual float32 keys; incomplete loop."""
import hashlib
import json
from pathlib import Path
import sys
from physical_length_window import ROOT,UNIT
sys.path.insert(0,str(UNIT/'native-qualification'))
import native_curve_gate as gate

def run():
    path=ROOT/'arrays/trial_968_prefix_native.npz'
    print(json.dumps({'stage':'native_prefix_gate_started','input':str(path)}),flush=True)
    report=gate.qualify_npz(path,identity_world=True,
                           config=gate.GateConfig(wall_seconds=60,max_leaves=250000,memory_limit_mib=900))
    report.update(scope='Native geometry of rejected loop968 body prefix and explicitly regenerated full-guide wrappers only.',
                  production_loop_complete=False,three_loop_qualified=False,selected_r5_promoted=False)
    output=ROOT/'receipts/native_prefix_gate.json'
    output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'qualified':report['qualified'],'status':report['status'],
                      'elapsed_seconds':report.get('elapsed_seconds'),'failures':report['failures'][:2]}),flush=True)

if __name__=='__main__':run()
