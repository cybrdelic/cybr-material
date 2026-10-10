# Approved Victorville ground R5D

Material 36 selects the ground appearance approved on 2026-10-09. The accepted construction is 73,733 separate mineral particles over a compacted fine substrate in a finite 0.18 m specimen. Three particle bands span 0.16–5.8 mm. This is Blender geometry plus PBR shading, not an equivalent flat map material or a qualified seamless tile.

The three files in `src/` are byte-identical to the delivered R5D checkpoint. `PROVENANCE.json` preserves the checkpoint's original work-in-progress status and separately records the later visual approval. R2 remains rejected research; later road scenes, distant LOD bakes and asphalt substitutions are not selected here. The existing approved asphalt selection is unchanged.

## Reconstruct in a fresh directory

Dependencies: Python with NumPy/SciPy and Blender 4.3.2. From this folder:

    python reconstruct.py --out /path/to/new/r5d-build --prepare-only
    python reconstruct.py --out /path/to/another/new/r5d-build --blender blender

The first command regenerates and verifies the fracture prototypes and microbed arrays. The second additionally runs the original scene builder and creates `preview/diagnostic_d.blend`. Neither command renders. They refuse an existing output directory and never write into this pinned source bundle. The scene builder uses two CPU threads. Do not overlap full geometry reconstruction with memory-heavy renders on a small host.

`expected_inputs.json` verifies generated array contents, not timestamp-bearing NPZ container bytes. `receipts/geometry_d.json` is the original execution receipt, not a new test result. `REVIEW.json` pins the original accepted scene and image hashes. Source-only CI tests reconstruction inputs and selection boundaries; it does not claim a fresh Blender rebuild, image equivalence, measured material calibration or full formation physics.

The registry's strict selected-scene build/capture path requires the original accepted scene at `preview/diagnostic_d.blend` with its pinned hash. A newly saved reconstruction may have different Blender binary bytes and does not silently replace that exact selection. Rendered images, large Blender scenes and generated NPZ inputs are intentionally absent from this source publication.

The original scene uses analytic daylight, Cycles CPU, 48 fixed samples and AgX Medium High Contrast. Guided review output used OpenImageDenoise 2.5.1. No image pixels, scans or image-generation outputs feed the material.
