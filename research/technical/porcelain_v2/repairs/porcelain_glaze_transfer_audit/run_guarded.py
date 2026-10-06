"""Run one Blender inspection, terminate above 600 MiB RSS or 120 seconds."""
import json, os, signal, subprocess, sys, time
from pathlib import Path

root = Path(__file__).resolve().parent
command = ['/usr/bin/blender', '-b', '-t', '1', '--factory-startup', '--python',
           str(root / 'inspect_scene.py'), '--', sys.argv[1], str(root / 'scene_data.json')]
start = time.monotonic()
peak_kib = 0
with (root / 'inspection.log').open('w') as log:
    process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
                               env={**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'})
    while process.poll() is None:
        try:
            status = Path(f'/proc/{process.pid}/status').read_text()
            rss = next((int(s.split()[1]) for s in status.splitlines() if s.startswith('VmRSS:')), 0)
        except FileNotFoundError:
            rss = 0
        peak_kib = max(peak_kib, rss)
        if rss > 600 * 1024 or time.monotonic() - start > 120 or (root / 'inspection.log').stat().st_size > 1024 * 1024:
            os.killpg(process.pid, signal.SIGKILL)
            process.wait()
            break
        time.sleep(.1)
receipt = {'command': command, 'returncode': process.returncode,
           'peak_rss_mib': peak_kib / 1024, 'elapsed_seconds': time.monotonic() - start,
           'rss_limit_mib': 600, 'wall_limit_seconds': 120}
(root / 'resources.json').write_text(json.dumps(receipt, indent=2))
print(json.dumps(receipt, indent=2))
sys.exit(process.returncode)
