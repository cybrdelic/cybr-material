# Reviewed asphalt fracture source · 2026-10-09

The explicitly approved asphalt-only integration selects material 28's reviewed fracture model. It retains the previous selection and does not select the subsequent binder refinement or change other material identities.

## Reproduce

Use Python 3.12 and the material-production requirements.txt. From this directory:

    OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python source/generate_candidate.py
    python source/verify_native.py
    python source/verify_preservation.py
    python source/test_height_normals.py

All four native4096 maps and GeometryHeight.npy are evaluated directly from procedural source. No maps, scans or packed scene are required as input. Generation uses approximately 4.1 GiB peak RAM; run serially. expected_native.json pins the exact reviewed output hashes.

The material-production command also supports a bounded smoke evaluation:

    python -m materials generate 28 --resolution 64 --out runs/asphalt_smoke_unique

Run that command from material-production. Reduced-resolution outputs do not replace the reviewed native assets. Then generate the reviewed coupon and studio from the included numeric recipe:

    blender -b -t 2 --python-exit-code 17 --python source/build_scene.py

Use Blender 4.3.2. This creates CYBR_Asphalt_Fracture_Native4096_Candidate.blend at the registry's relative path. No external scene is an input. The builder and unified build inspection compare a render-relevant identity covering geometry, UVs, material bindings/nodes/map bytes, camera, lights, world and display. Saved blend-file bytes can vary with paths; the registry keeps the original reviewed binary hash as provenance while selected-build acceptance uses the independently inspected render contract. A changed scene cannot pass merely by editing a sidecar receipt.

## Representation and limits

The source preserves coarse seeded placements and finite irregular outlines, replaces angular domes with finite-bevel faceted crown precursors under nonlinear compaction, and adds fine mineral exposure coupled to geometry. Full metric OBJECT-space RGB16 normals replace residual tangent transfer, with strength 1 and Non-Color sampling. Geometry uses a 512-cell macro grid over the 0.25 m tile; the full normal derives from the complete height field using the 0.004 m height scale.

OBJECT normals support tested rotations, nonuniform scales and linked instances on this static coupon. Deformation or chart changes require regeneration. Fine normal relief does not fully represent self-shadowing geometry. Direct shader-output proof was before the BSDF grazing-normal correction.

This is an authored approximation, not mechanically solved compaction or measured bitumen/mineral optics. Bright uniform beadlike fines and polygonal plate-like crowns remain. Review was a matched native detail comparison; it does not establish universal realism or a new complete-panel selection for other materials.
