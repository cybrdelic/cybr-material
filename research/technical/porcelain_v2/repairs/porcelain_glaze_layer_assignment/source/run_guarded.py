"""Bound the admitted source-only Blender process. Never call a renderer."""
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

root = Path(__file__).resolve().parents[1]
out = root / 'scenes'
out.mkdir(exist_ok=True)
command = ['/usr/bin/blender', '-b', '-t', '1', '--factory-startup', '--python-exit-code', '2',
           '--python', str(root / 'source/build_layer_candidate.py'), '--', '--admitted']
start, peak, reason = time.monotonic(), 0, None
with (out / 'build.log').open('w') as log:
    process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
                               env={**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'})
    while process.poll() is None:
        try:
            status = Path(f'/proc/{process.pid}/status').read_text()
            rss = next((int(s.split()[1]) for s in status.splitlines() if s.startswith('VmRSS:')), 0)
        except FileNotFoundError:
            rss = 0
        peak = max(peak, rss)
        size = sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
        if rss > 1024 * 1024 or time.monotonic() - start > 120 or size > 8 * 1024 * 1024:
            reason = 'Resource ceiling exceeded'
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            break
        time.sleep(.05)
receipt = {'command': command, 'returncode': process.returncode, 'peak_rss_mib': peak / 1024,
           'elapsed_seconds': time.monotonic() - start, 'rss_limit_mib': 1024,
           'wall_limit_seconds': 120, 'output_limit_mib': 8, 'stop_reason': reason,
           'rendered': False}
(out / 'resources.json').write_text(json.dumps(receipt, indent=2))
print(json.dumps(receipt, indent=2))
sys.exit(process.returncode)
