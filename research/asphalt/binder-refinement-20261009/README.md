# CYBR Asphalt · binder refinement, compact source

This is the unselected 2026-10-09 graded-fine/binder-exposure candidate. It preserves the reviewed fractured coarse-grain construction and replaces the uniform exposed fine-mineral field with three size bands, geometry-derived channel grading, binder-contact burial and partial exposed-crown optics.

This source-only package regenerates the native maps, packed scene and exact matched capture job. No image, scan, input `.blend`, private artifact or external geometry is needed. `scene_recipe.json` contains the numerical studio, shader and mesh construction. `expected_native.json` pins the generated image/geometry bytes; `expected_scene.json` pins the actual candidate's render-relevant scene signature.

## Reproduce

Use Python 3.12, requirements.txt, Blender 4.3.2, and official Intel OIDN 2.5.1. Run serially from the extracted package directory. Set `OPENBLAS_NUM_THREADS=2`, `OMP_NUM_THREADS=2`, `MKL_NUM_THREADS=2`, and `NUMEXPR_NUM_THREADS=2`. Synthesis needs approximately 4 GiB RAM; the matched beauty capture has a 3 GiB cap.

1. `python3 source/test_refinement.py`
2. `python3 source/generate_candidate.py`
3. `python3 source/verify_native.py`
4. Optional numerical normal check: `python3 source/verify_normal_encoding.py`
5. `blender -b -t 2 --python-exit-code 17 --python source/build_scene.py`
6. Set `OIDN_BIN` to the official `oidnDenoise` executable, then run `bash pipeline/clean_render.sh render/job.json`
7. Optional actual shader-normal proof: `blender -b -t 2 --python-exit-code 17 --python source/verify_renderer_normals.py -- native`

The builder writes `CYBR_Asphalt_Binder_Refinement_Native4096_Candidate.blend` and regenerates `render/job.json` with current paths and the packed scene's hash. The clean matched image is `render/CLEAN_OIDN_PANEL.png`.

## Mechanism and units

Native 4096×4096, tile 0.25 m, height `(encoded - 0.5) * 0.004 m`, 513×513 macro grid. The RGB16 normal map is a full metric OBJECT-space normal for this flat coupon, not a tangent residual or a general curved-mesh adapter.

The fine bands have nominal apothems of 0.085–0.170 mm, 0.140–0.290 mm, and 0.240–0.465 mm. Coarse geometry determines contact support and channel width, which influence fine-size grading, binder menisci and burial. The crown remaining above its binder contact drives relief, exposed-mineral color and roughness with one surface owner, including antialiased edges. No independent matrix color/roughness noise is added. Reviewed coarse microrelief is retained and is not optically coupled to matrix exposure underneath a fully covered coarse grain.

## Validation and visual limits

The compact builder was independently run against the existing candidate maps, without resynthesizing them. Its render signature exactly matched the already rendered candidate's geometry, UVs, materials, bound images, camera, lights, world and display. The generated capture profile matched the actual 960×960, 128-sample Cycles CPU/OIDN render. The source ZIP intentionally excludes all map, scene, pass and preview binaries.

The original candidate's native-integrity, coarse-preservation and normal checks passed. The actual shader-normal AOV maximum vector error was 1.46e-7 across 102,400 pixels. The visual improvement is reduced uniform fine-bead sparkle and more continuous binder; some smooth binder pockets and the preserved coarse plate-like appearance remain. This candidate is unselected and not physically qualified. It is authored procedural appearance, not measured asphalt, solved particle packing, force-balanced compaction or calibrated binder flow. Fine self-occlusion is limited by the normal-map representation. The compact scene was signature-verified rather than rerendered.
