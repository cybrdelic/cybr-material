"""One-core, CPU-only subprocess watchdog for scoped numerical jobs."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import psutil

def run():
    parser=argparse.ArgumentParser()
    parser.add_argument('--seconds',type=float,default=90)
    parser.add_argument('--memory-mib',type=float,default=1000)
    parser.add_argument('--receipt',type=Path,required=True)
    parser.add_argument('command',nargs=argparse.REMAINDER)
    args=parser.parse_args()
    command=args.command[1:] if args.command[:1]==['--'] else args.command
    if not command:parser.error('missing command')
    env=os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',
               NUMEXPR_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',CUDA_VISIBLE_DEVICES='',
               PYTHONPATH=str(Path(__file__).resolve().parent))
    available=psutil.virtual_memory().available
    if available<2*1024**3:raise RuntimeError('Less than2GiB available; numerical job not started')
    start=time.monotonic()
    log=args.receipt.with_suffix('.log')
    log.parent.mkdir(parents=True,exist_ok=True)
    peak=0;status='completed';cpu_seconds=0
    with log.open('w',encoding='utf8') as stream:
        child=subprocess.Popen(command,env=env,stdout=stream,stderr=subprocess.STDOUT)
        process=psutil.Process(child.pid)
        affinity=process.cpu_affinity() if hasattr(process,'cpu_affinity') else []
        if affinity:process.cpu_affinity([affinity[-1]])
        if os.name=='nt':process.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
        while child.poll() is None:
            try:
                processes=[process]+process.children(recursive=True)
                rss=sum(p.memory_info().rss for p in processes if p.is_running())
                peak=max(peak,rss)
                times=[p.cpu_times() for p in processes if p.is_running()]
                cpu_seconds=max(cpu_seconds,sum(x.user+x.system for x in times))
                reason=('wall_budget_exceeded' if time.monotonic()-start>args.seconds else
                        'memory_budget_exceeded' if rss>args.memory_mib*1024**2 else
                        'system_memory_pressure' if psutil.virtual_memory().available<1.5*1024**3 else None)
                if reason:
                    status=reason
                    for p in reversed(processes):
                        try:p.kill()
                        except psutil.NoSuchProcess:pass
                    break
            except psutil.NoSuchProcess:pass
            time.sleep(.1)
        code=child.wait(timeout=5)
    result={'status':status,'exit_code':code,'elapsed_seconds':time.monotonic()-start,
            'peak_process_tree_rss_mib':peak/1024**2,'observed_process_tree_cpu_seconds':cpu_seconds,
            'initial_available_memory_mib':available/1024**2,'cpu_affinity':[affinity[-1]] if affinity else None,
            'numerical_threads':1,'gpu_disabled':True,'memory_limit_mib':args.memory_mib,
            'wall_limit_seconds':args.seconds,'command':command,'log':str(log)}
    args.receipt.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
    return code if status=='completed' else 3

if __name__=='__main__':raise SystemExit(run())
