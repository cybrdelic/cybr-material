"""One admitted numerical probe, bounded to 1 GiB RSS and 180 seconds."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time,resource

ROOT=Path(__file__).resolve().parents[1]
receipt=ROOT/'receipts/fixed_budget_3d_probe.json'
if receipt.exists():
    raise SystemExit('Existing 3D receipt must be inspected; refusing duplicate run.')
(ROOT/'arrays').mkdir(exist_ok=True)
env=os.environ.copy()
env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
start=time.monotonic();peak=0;reason=None
source={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'source').glob('*.py')}
with (ROOT/'receipts/fixed_budget_3d_probe.log').open('x') as log:
    proc=subprocess.Popen([sys.executable,'-u',str(ROOT/'source/probe_3d_correction.py')],env=env,stdout=log,stderr=subprocess.STDOUT)
    while proc.poll() is None:
        try:
            status=Path(f'/proc/{proc.pid}/status').read_text()
            rss=int(next(line.split()[1] for line in status.splitlines() if line.startswith('VmRSS:')))*1024
            peak=max(peak,rss)
        except (FileNotFoundError,StopIteration):pass
        if peak>1024**3:reason='RSS cap exceeded'
        if time.monotonic()-start>180:reason='Wall cap exceeded'
        if reason:
            proc.kill();proc.wait();break
        time.sleep(.1)
result={'exit_code':proc.returncode,'guard_stop':reason,'wall_seconds':time.monotonic()-start,'sampled_peak_rss_mib':peak/1024**2,'kernel_child_peak_rss_mib':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/1024,'rss_cap_mib':1024,'wall_cap_seconds':180,'source_sha256':source,'receipt_exists':receipt.exists()}
(ROOT/'receipts/fixed_budget_3d_execution.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result))
raise SystemExit(proc.returncode if proc.returncode is not None else 1)
