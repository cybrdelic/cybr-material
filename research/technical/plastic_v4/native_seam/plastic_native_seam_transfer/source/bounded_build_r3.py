from pathlib import Path
import subprocess,os,signal,time,json,resource
R=Path(__file__).resolve().parents[1];start=time.monotonic();peak=0;reason=None
with (R/'receipts/build_r3.log').open('x') as f:
 p=subprocess.Popen(['/usr/bin/blender','-b','-t','1','--python-exit-code','1','--python',str(R/'source/build_fixture_r3.py')],stdout=f,stderr=subprocess.STDOUT,start_new_session=True,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
 while p.poll() is None:
  try:status=Path(f'/proc/{p.pid}/status').read_text();rss=int(next(x for x in status.splitlines() if x.startswith('VmRSS:')).split()[1])/1024
  except (OSError,StopIteration):rss=0
  peak=max(peak,rss)
  if rss>600:reason='RSS cap exceeded'
  if time.monotonic()-start>180:reason='Build time cap exceeded'
  if reason:
   os.killpg(p.pid,signal.SIGTERM)
   try:p.wait(5)
   except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL)
   break
  time.sleep(.1)
 code=p.wait()
r={'exit_code':code,'sampled_peak_MiB':peak,'kernel_peak_MiB':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/1024,'seconds':time.monotonic()-start,'stop_reason':reason,'cap_MiB':600,'scope':'One source build; no render'};(R/'receipts/build_performance_r3.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r));raise SystemExit(code)
