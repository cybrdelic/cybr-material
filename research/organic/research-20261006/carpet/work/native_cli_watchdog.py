"""Lightweight Windows CLI wall watchdog, including worker import time."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time

def main(script,arguments):
    parser=argparse.ArgumentParser(add_help=False)
    parser.add_argument('--wall-seconds',type=float,default=180)
    parser.add_argument('--output',type=Path,required=True)
    args,_=parser.parse_known_args(arguments)
    start=time.monotonic()
    child=subprocess.Popen([sys.executable,str(script),*arguments,'--_worker'])
    while child.poll() is None:
        if time.monotonic()-start>=max(.1,args.wall_seconds):
            child.kill()
            child.wait(timeout=5)
            result={'qualified':False,'complete':False,
                    'status':'unproven_cli_hard_wall_budget_exceeded',
                    'failures':[{'code':'cli_hard_wall_budget_exceeded'}],
                    'elapsed_seconds':time.monotonic()-start,
                    'wall_includes_worker_imports':True}
            args.output.parent.mkdir(parents=True,exist_ok=True)
            args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
            print(json.dumps({'qualified':False,'status':result['status'],'receipt':str(args.output)}))
            return 2
        time.sleep(.05)
    return child.returncode
