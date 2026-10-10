# CYBR warm grazing preview preset

`cybr_grazing_warm_v1` is the local default for **new material previews prepared with this entry point**. It captures the user's chosen corrected-asphalt grazing presentation. This is an additive presentation adapter, not a change to an accepted material or its physical construction.

The exact numeric values are in `cybr_grazing_warm_v1.json`; `preview_defaults.json` selects it. Historical capture jobs and the existing render pipeline remain unchanged. No old jobs or the full material collection are automatically rerendered. No Blender user preferences or operating-system settings are saved.

## Prepare a preview without rendering

Requires Blender 4.3.2 (the inspected reference version) and Python 3. All output goes into a **new, nonexistent directory whose parent already exists**.

```sh
python prepare_preview.py \
  --source /path/to/immutable_source.blend \
  --out /path/to/new_preview_directory \
  --specimen-width-m 0.25 \
  --specimen-pivot-m 0 0 0.008
```

The default is `--composition reference`. Outputs are `preview.blend`, `job.json`, and `preparation_receipt.json`. This command **cannot render**. It launches Blender with automatic embedded scripts disabled, loads the source read-only in practice, applies presentation settings in memory, and saves only the derivative. It hashes the source before/after and refuses existing output directories, symlinks, or overwriting the source. Optional `--expect-source-sha256 SHA` pins the input. Repeat `--protect-file PATH` to verify registries or other important files remain unchanged. Failures may leave a partial fresh output directory; inspect it and choose a new output directory rather than overwriting it.

For repo integration, this directory can live at `material-production/preview_presets/`; run `python material-production/preview_presets/prepare_preview.py ...`. Existing source-faithful audit/capture commands retain their existing meaning. Future presentation previews should use this entry point explicitly.

## Physical-scale contract

- The specimen's actual geometry, transforms, physical feature sizes, material assignments, shader graphs, and textures are not rescaled or rebuilt.
- The source must use metre coordinates with Blender unit scale 1. Other unit scales fail closed. Dimensions must be supplied from the specimen, not inferred from a filename or a surrounding floor.
- `specimen_width_m` is the **actual physical specimen width** corresponding to the reference tile's 0.25 m width.
- `specimen_pivot_m` is the **world-space surface-centre datum** corresponding to `(0, 0, 0.008)` in the asphalt reference, not necessarily an object's origin or bounding-box centre. For a surface centred on world z=0, supply `0 0 0`.
- For scale `s = specimen_width_m / 0.25`, rig locations transform as `new_pivot + (reference_location - reference_pivot) * s`. Emitter diameter and camera distances/clips scale by `s`, key power by `s²`. Source geometry is untouched. World intensity and light colour remain unchanged.
- For a 0.60 m specimen this gives `s=2.4`, key diameter about 0.132 m and power about 18.6624 W. Use the JSON values rather than rounding these examples.
- This is a Z-up, reference-aligned rig. A specimen in a different orientation requires an explicit future adapter; the helper never silently rotates geometry.

## Composition modes

**Reference (default):** reproduces the original relative camera position, recorded XYZ Euler rotation, orthographic span, shifts, and crop. At width 0.25 m and pivot `(0,0,0.008)`, it uses the exact reference camera values. Crucially, the reference camera does **not** point exactly at the light's target, and its orientation is not recomputed. This is a detail crop of the original tile, not a promise that the full specimen is visible.

**Full specimen:** keeps the same light and camera orientation but moves the camera in its image plane and fits a padded view around only the explicitly named specimen objects. Name every part that should be included; do not name the backdrop or floor.

```sh
python prepare_preview.py --source /path/source.blend --out /path/new_full \
  --specimen-width-m 0.60 --specimen-pivot-m 0 0 0 \
  --composition full --specimen-object 'Ground mesh' \
  --specimen-object 'Gravel mesh' --padding 1.10
```

Alternatively supply `--specimen-bounds-m MIN_X MIN_Y MIN_Z MAX_X MAX_Y MAX_Z`, six explicit world-space metre coordinates. Do not combine bounds with object names. Bounding-box fits are conservative; procedural modifiers whose evaluated bounds differ need explicit evaluated bounds from the caller.

**Detail:** changes framing only, with an explicit horizontal view width. Crop centre defaults to the surface pivot. Camera orientation and lights are preserved.

```sh
python prepare_preview.py --source /path/source.blend --out /path/new_detail \
  --specimen-width-m 0.60 --specimen-pivot-m 0 0 0 \
  --composition detail --crop-width-m 0.12 --crop-center-m 0.03 -0.02 0
```

`--width` and `--height` specify output pixels (default 640×640); rectangular images are supported. Full/detail composition modifies framing only, not the light, geometry or feature scale. Detail view width is not a specimen-size override.

## Blender API

```python
from preview_presets import apply_preset, build_job
plan = apply_preset(
    bpy.context.scene,
    specimen_width_m=0.60,
    specimen_pivot_m=(0, 0, 0),
    composition='reference',
)
```

`apply_preset` only changes the loaded scene in memory; it never saves or renders. The source-preserving CLI should be preferred whenever a derivative file is needed. `plan_preset` can inspect numeric settings without importing Blender. Full composition accepts `specimen_objects=['exact object name', ...]` through `apply_preset`, or `specimen_bounds_m=[[min_x,min_y,min_z],[max_x,max_y,max_z]]` through either function. Detail accepts `crop_width_m` and optional `crop_center_m`.

## Exact capture contract

One warm DISK area light at the recorded grazing position; other explicit scene lights have energy zero. Constant world approximately RGBA `(0.18,0.21,0.25,1)`, strength `0.04`. AgX / Medium High Contrast, exposure 0, gamma 1. Cycles CPU, 128 samples, seed 1810, fixed sampling, no light tree, no native denoising; bounces and clamps copied from the receipt. Render output is native RGB16 PNG, dither zero. The generated job requests linear HDR beauty plus native albedo/normal guides for Intel Open Image Denoise 2.5.1 CPU RT high-quality guided denoising.

To actually render a prepared job, use the separately approved existing pipeline, with its OIDN 2.5.1 executable available:

```sh
OIDN_BIN=/path/to/oidn-2.5.1/bin/oidnDenoise \
  /path/to/pipeline/clean_render.sh /path/to/new_preview_directory/job.json
```

Preparing the scene does not verify that OIDN is installed or produce a render. The job records the derivative source hash and preset hash. The existing pipeline produces its own verified pipeline snapshot and final capture receipt. Compare that receipt's pipeline hash to the recorded reference pipeline hash when exact pipeline equivalence matters. There is no colour-only denoising fallback in this contract. Daylight is an optional separately requested secondary presentation, never the default here.

## Review limitations

- **This presentation changes lighting.** Do not call an old-versus-new image comparison “matched” unless both use the same preset version, physical rig scale, world, camera/framing, exposure, renderer and capture pipeline.
- Matching rig settings alone does not prove identical pixels, material quality, or visual acceptance. Physical specimen scale must be exact and separately verified.
- It is a presentation default, not a claim that warm grazing light is the best neutral colour/optical diagnostic for every material.
- Emissive source materials and geometry remain untouched. Such sources may emit additional light. A single explicit light is not a guarantee of a single illumination contribution.
- Source compositing/sequencer effects and render borders are disabled in the derivative so they do not contaminate the new preview. Material shaders and texture content are preserved. New camera and world datablocks isolate the rig from source camera constraints or HDRI worlds.
- Static frame 1 previews are supported. Animated/deformed sources require separately qualified evaluation; the CLI rejects detected geometry-transform/material-slot changes while preparing.
- Tested on Blender 4.3.2. Missing required properties fail rather than silently relaxing the preset. Different Blender/OIDN versions require renewed comparison.
- These checks perform no rendering. No visual-equivalence or cross-material lighting claim is made.

## Verification

```sh
python -m unittest discover -s /path/to/preview_presets -v
```

Tests cover default identity, exact reference camera (including its non-target aim), uniform scale/power scaling, full/detail framing, invalid inputs, new-output/source-hash guards, and the no-render entry point. A Blender source-to-derivative preparation additionally checks recorded settings and immutable-source/registry hashes without rendering.

The CLI reopens the saved derivative before recording its inspected settings. To verify its settings and immutable-file hashes afterward:

```sh
python verify_preparation.py /path/to/new_preview_directory/preparation_receipt.json
```
