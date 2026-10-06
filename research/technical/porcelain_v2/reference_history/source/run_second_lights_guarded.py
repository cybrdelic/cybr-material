"""One admitted second-rig derivative operation, 1 GiB /120 s /4 MiB."""
import json,os,signal,subprocess,sys,time
from pathlib import Path
R=Path(__file__).resolve().parents[1]; O=R/'second_light_sources';O.mkdir(exist_ok=True)
assert '--admitted' in sys.argv
receipt=O/'resources.json'; assert not receipt.exists()
command=['/usr/bin/blender','-b','-t','1','--factory-startup','--python-exit-code','2','--python',str(R/'source/build_second_lights.py'),'--','--admitted']
start=time.monotonic();peak=0;reason=None
with (O/'build.log').open('x') as log:
    proc=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'})
    while proc.poll() is None:
        try:rss=next((int(s.split()[1]) for s in Path(f'/proc/{proc.pid}/status').read_text().splitlines() if s.startswith('VmRSS:')),0)
        except FileNotFoundError:rss=0
        peak=max(peak,rss);size=sum(p.stat().st_size for p in O.rglob('*') if p.is_file());elapsed=time.monotonic()-start
        if rss>1024*1024 or elapsed>120 or size>4*1024**2:
            reason={'rss_kib':rss,'elapsed_s':elapsed,'output_bytes':size};os.killpg(proc.pid,signal.SIGKILL);proc.wait();break
        time.sleep(.05)
result={'command':command,'returncode':proc.returncode,'peak_observed_rss_mib':peak/1024,'elapsed_seconds':time.monotonic()-start,'output_bytes':sum(p.stat().st_size for p in O.rglob('*') if p.is_file()),'rss_limit_mib':1024,'wall_limit_seconds':120,'output_limit_mib':4,'stop_reason':reason,'rendered':False}
receipt.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));sys.exit(proc.returncode)
