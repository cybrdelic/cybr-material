"""Explicit isolated build. No installation and no production library overwrite."""
from pathlib import Path
import hashlib, json, platform, resource, subprocess, time
ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'source/directional_native.cpp'
LIBRARY = ROOT/'source/libdirectional_native.so'
FLAGS = ['-std=c++17','-O3','-fPIC','-shared','-fno-fast-math','-ffp-contract=off','-Wall','-Wextra','-Wpedantic']
def main():
    start=time.perf_counter()
    compiler=subprocess.check_output(['g++','--version'],text=True).splitlines()[0]
    command=['g++',*FLAGS,str(SOURCE),'-o',str(LIBRARY)]
    subprocess.run(command,check=True,timeout=100)
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    receipt={'compiler':compiler,'flags':FLAGS,'platform':platform.platform(),'source_sha256':sha(SOURCE),'library_sha256':sha(LIBRARY),'seconds':time.perf_counter()-start,'maximum_child_resident_memory_KiB':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss}
    (ROOT/'receipts/build.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
