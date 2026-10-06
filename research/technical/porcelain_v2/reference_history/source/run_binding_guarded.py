"""Admitted source binding only, bounded to 1 GiB /120 s /4 MiB new outputs."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'scenes'
assert '--admitted' in sys.argv
OUT.mkdir(exist_ok=True)
receipt=OUT/'resources.json'
assert not receipt.exists()
command=['/usr/bin/blender','-b','-t','1','--factory-startup','--python-exit-code','2',
         '--python',str(ROOT/'source/bind_sources.py'),'--','--admitted']
start=time.monotonic()
peak=0
reason=None
with (OUT/'build.log').open('x') as log:
    proc=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,
             env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'})
    while proc.poll() is None:
        try:
            status=Path(f'/proc/{proc.pid}/status').read_text()
            rss=next((int(s.split()[1]) for s in status.splitlines() if s.startswith('VmRSS:')),0)
        except FileNotFoundError: rss=0
        peak=max(peak,rss)
        size=sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file())
        elapsed=time.monotonic()-start
        if rss>1024*1024 or elapsed>120 or size>4*1024**2:
            reason={'rss_kib':rss,'elapsed_s':elapsed,'output_bytes':size}
            os.killpg(proc.pid,signal.SIGKILL)
            proc.wait()
            break
        time.sleep(.05)
result={'command':command,'returncode':proc.returncode,'peak_observed_rss_mib':peak/1024,
        'elapsed_seconds':time.monotonic()-start,'output_bytes':sum(p.stat().st_size for p in OUT.rglob('*') if p.is_file()),
        'rss_limit_mib':1024,'wall_limit_seconds':120,'output_limit_mib':4,'stop_reason':reason,
        'rendered':False,'selected':False}
receipt.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
sys.exit(proc.returncode)
