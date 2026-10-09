"""Bounded single build: admission must already be granted by the parent."""
from pathlib import Path
import subprocess,os,signal,time,json,resource,sys
R=Path(__file__).resolve().parents[1];name,cap,seconds=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]);assert name in {'build_shader_fixture.py','build_complete_candidate.py'};script=R/'source'/name;tag=name.removesuffix('.py')+('_'+sys.argv[4] if len(sys.argv)>4 else '');start=time.monotonic();peak=0;reason=None
with (R/'receipts'/f'{tag}.log').open('x') as f:
    p=subprocess.Popen(['/usr/bin/blender','-b','-t','1','--python-exit-code','1','--python',str(script)],stdout=f,stderr=subprocess.STDOUT,start_new_session=True,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
    while p.poll() is None:
        try:status=Path(f'/proc/{p.pid}/status').read_text();rss=int(next(x for x in status.splitlines() if x.startswith('VmRSS:')).split()[1])/1024
        except (OSError,StopIteration):rss=0
        peak=max(peak,rss)
        if rss>cap:reason='RSS cap exceeded'
        if time.monotonic()-start>seconds:reason='Build time cap exceeded'
        if reason:
            os.killpg(p.pid,signal.SIGTERM)
            try:p.wait(5)
            except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL)
            break
        time.sleep(.1)
    code=p.wait()
r={'exit_code':code,'sampled_peak_MiB':peak,'kernel_peak_MiB':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/1024,'seconds':time.monotonic()-start,'stop_reason':reason,'cap_MiB':cap,'scope':'One admitted source build; no render'};(R/'receipts'/f'{tag}_performance.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r));raise SystemExit(code)
