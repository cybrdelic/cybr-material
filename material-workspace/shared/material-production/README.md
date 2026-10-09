# Selected material pipeline

The registry contains 35 identities with exact source/capture hashes. It preserves
stronger visual anchors and never selects the newest experiment automatically.
The runtime and pinned original/aggregate generators are included; native scenes,
map sets and some per-material reconstruction inputs require exact restoration.

## Source-only use

From the repository root:

```sh
python -m pip install -r material-workspace/shared/material-production/requirements.txt
python audit/verify_checkpoint.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python audit/run_portable_tests.py
cd material-workspace/shared/material-production
python -m materials list
python -m materials doctor
python -m materials generate 05 --resolution 64 --out runs/brass_fields_unique
python -m materials generate 32 --resolution 64 --out runs/concrete_fields_unique
```

`doctor` distinguishes missing assets from quality acceptance. The source-only test
runner explicitly excludes one leather native-scene check; it does not turn that
missing asset into a pass. Reduced-resolution generation is a numerical smoke test,
not native 4K export or appearance evidence. Each output directory must be new.

## Build and capture after restoring assets

Restore the exact scene/map versions and hashes referenced by `registry.json` and
[`audit/ARTIFACTS.json`](../../../audit/ARTIFACTS.json). Preserve the
`material-workspace/shared` directory layout. Set `MATERIAL_WORKSPACE` to the absolute
`material-workspace` directory and use `MATERIAL_PATH_MAP_FILE` for legacy embedded
paths. These settings do not create missing assets or certify full portability.

Captures require Blender 4.3.2 and OIDN 2.5.1 in their configured locations. The shared
capture implementation is in `../cloud-eevee`. After restoring dependencies:

```sh
python -m materials build 05 --out runs/brass_build_unique
python -m materials capture runs/brass_build_unique/build.json --out runs/brass_capture_unique
```

`build` verifies selected scene/map hashes and inspects real geometry, UVs, scale and
shader bindings. `capture` revalidates dependencies and runs the shared Cycles → EXR
→ selected filtering mode → RGB16 pipeline in a new output directory. It does not
overwrite a selected scene or promote a candidate. See [capture modes](CAPTURE_MODES.md).

## Qualification boundary

- The brass and polished-concrete generation/inspection paths have focused
  integration coverage. Other catalogue entries do not inherit that qualification.
- Native brass decoded pixels match, but historical compressed PNG bytes differ.
  Concrete native channel bytes match; concrete appearance remains unresolved.
  Exact receipts remain necessary for strict reconstruction/parity workflows.
- Wood anchors and r5 fleece/carpet remain protected. New numerical or structural
  studies do not replace a stronger reviewed complete panel.
- Full formation physics, measured calibration and clean-clone reconstruction of
  all 35 materials remain incomplete. Historical iterations are preserved separately.


## Selected asphalt update · 2026-10-09

Material28 now selects the explicitly reviewed fracture candidate. [Exact source, native generation and restoration](vendor/asphalt-reviewed-20261009/README.md) retain the previous selection separately. This asphalt-only update does not select the subsequent binder refinement or change other material identities.
