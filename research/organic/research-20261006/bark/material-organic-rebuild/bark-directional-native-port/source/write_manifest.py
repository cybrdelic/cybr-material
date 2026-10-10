"""Record exact port, frozen-source and loaded-runtime dependency hashes."""
from pathlib import Path
import hashlib,importlib.metadata,json,platform,shutil,subprocess,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    sources={str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'source').iterdir()) if p.is_file()}
    references=[BASE/'bark-constitutive-upgrade/source/directional_candidate.py',BASE/'bark-constitutive-upgrade/source/quadratic_prism_frozen.py',BASE/'bark-constitutive-incremental/source/incremental_envelope.py',BASE/'causal-bark/volume-core/source/batch_volume_incremental.cpp',BASE/'causal-bark/volume-core/source/native_volume.py']
    dependencies={str(p.relative_to(BASE)):sha(p) for p in references}
    runtime_paths={Path(sys.executable).resolve(),Path(shutil.which('g++')).resolve(),Path(np.__file__)}
    numpy_root=Path(np.__file__).parent
    runtime_paths.update(numpy_root.rglob('*.so'))
    library=ROOT/'source/libdirectional_native.so'
    ldd=subprocess.check_output(['ldd',str(library)],text=True)
    for line in ldd.splitlines():
        for token in line.split():
            if token.startswith('/') and Path(token).is_file():runtime_paths.add(Path(token).resolve())
    runtime={str(p):sha(p) for p in sorted(runtime_paths)}
    receipts={str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'receipts').iterdir()) if p.is_file() and p.name!='manifest.json'}
    result=dict(status='ISOLATED_SYNTHETIC_PORT_ONLY',python=platform.python_version(),numpy=importlib.metadata.version('numpy'),platform=platform.platform(),documentation_sha256={'README.md':sha(ROOT/'README.md')},sources_sha256=sources,frozen_dependency_sha256=dependencies,runtime_binary_and_numpy_extension_sha256=runtime,native_link_dependencies=ldd,receipts_sha256=receipts,notes='Hashes identify the inspected source files, native output and interpreter/compiler/runtime binaries, including NumPy extension modules. This is not a wheel lockfile or a hash of all operating-system files. No SciPy runtime is required by this port.')
    (ROOT/'receipts/manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'source_files_hashed':len(sources),'frozen_sources_hashed':len(dependencies),'runtime_files_hashed':len(runtime),'receipt_files_hashed':len(receipts)},indent=2))
if __name__=='__main__':main()
