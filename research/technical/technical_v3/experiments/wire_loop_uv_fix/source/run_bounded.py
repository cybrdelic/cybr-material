"""Source-only build monitor; separate PID tree, measured RSS and exit receipt."""
import sys,subprocess,time,json,os,signal
from pathlib import Path
import psutil
cap=int(sys.argv[1])*1024**2; receipt=Path(sys.argv[2]); command=sys.argv[3:]
t=time.monotonic();p=subprocess.Popen(command,start_new_session=True);peak=0;exceeded=False
while p.poll() is None:
 try:
  root=psutil.Process(p.pid);rss=sum(q.memory_info().rss for q in [root,*root.children(recursive=True)] if q.is_running());peak=max(peak,rss)
  if rss>cap:
   exceeded=True;os.killpg(p.pid,signal.SIGTERM);break
 except psutil.NoSuchProcess:pass
 time.sleep(.1)
try: code=p.wait(timeout=10)
except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);code=p.wait()
receipt.write_text(json.dumps(dict(command=command,return_code=code,peak_rss_mib=peak/1024**2,cap_mib=cap/1024**2,cap_exceeded=exceeded,elapsed_s=time.monotonic()-t),indent=2))
raise SystemExit(97 if exceeded else code)
