"""Admitted 1 GiB / 120 second builder, 8 MiB output ceiling, no render."""
import json, os, signal, subprocess, sys, time
from pathlib import Path
root=Path(__file__).resolve().parent
source,out=map(lambda x:Path(x).resolve(),sys.argv[1:3])
out.mkdir(parents=True,exist_ok=True)
command=['/usr/bin/blender','-b','-t','1','--factory-startup','--python-exit-code','2',
         '--python',str(root/'build_control_pair.py'),'--',str(source),str(out),'--admitted']
start=time.monotonic();peak=0;reason=None
with (out/'build.log').open('w') as log:
    process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,
        env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'})
    while process.poll() is None:
        try:
            status=Path(f'/proc/{process.pid}/status').read_text()
            rss=next((int(s.split()[1]) for s in status.splitlines() if s.startswith('VmRSS:')),0)
        except FileNotFoundError:
            rss=0
        peak=max(peak,rss)
        total=sum(p.stat().st_size for p in out.iterdir() if p.is_file())
        if rss>1024*1024 or time.monotonic()-start>120 or total>8*1024*1024:
            reason='Resource ceiling exceeded';os.killpg(process.pid,signal.SIGKILL);process.wait();break
        time.sleep(.1)
receipt={'command':command,'returncode':process.returncode,'peak_rss_mib':peak/1024,
         'elapsed_seconds':time.monotonic()-start,'rss_limit_mib':1024,'wall_limit_seconds':120,
         'output_limit_mib':8,'stop_reason':reason}
(out/'resources.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))
sys.exit(process.returncode)
