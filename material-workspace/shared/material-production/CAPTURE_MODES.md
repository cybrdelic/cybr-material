# Capture modes

Use the same immutable experiment/build contract for every mode. A capture never selects or promotes an experiment. Use a new output directory each time.

## Guided Cycles (default)

```sh
python -m materials capture runs/BUILD/build.json --out runs/NEW_CAPTURE
```

The renderer retains linear HDR beauty, albedo and normal EXRs, checks the declared guide contract, and runs OIDN. Output is `CLEAN_OIDN_PANEL.png`. A bounded guide alone does not establish detail preservation: inspect the actual image and compare important native texture detail before accepting it.

## Explicit unfiltered Cycles

```sh
python -m materials capture runs/BUILD/build.json --out runs/NEW_CAPTURE --mode unfiltered --samples 512 --fixed-samples
```

This disables adaptive sampling and captures exactly the requested maximum sample count. It retains HDR beauty and auxiliary EXRs, including diagnostic records of invalid auxiliary ranges, but never invokes OIDN. The native 16-bit PNG must round-trip through the HDR colour bridge within one code value. Output is `UNFILTERED_PANEL.png`; visible-noise acceptance remains a separate review. A failed guided job never silently falls back to this mode.

## Preserve display transmission and emission

```sh
python -m materials capture runs/BUILD/build.json --out runs/NEW_CAPTURE --mode component-recombined --samples 512 --fixed-samples
```

This explicit mode retains transmission plus emission and a float64 accumulation residual. OIDN filters only the nonnegative diffuse, glossy, environment and volume remainder. The source must supply valid declared albedo/normal guides. A full physical-pass export, a small numerical reconstruction guard, exact within-capture float32 identity, immutable inputs, and the native16 colour bridge are required. Output is `COMPONENT_RECOMBINED_PANEL.png`.

Review the retained signal, filtered remainder and final image separately. Preserving a noisy transmitted signal does not make it clean. This route has been reviewed for the current LCD/CRT recipes; other recipes need their own signal/noise review. Source optics and emitted structure are not tuned by the mode.

With no CLI mode override, capture uses an explicitly declared profile mode; profiles without one remain guided. A failure never changes processing mode.

## Workbench construction views

Register `engine: BLENDER_WORKBENCH` in the immutable experiment and use `--mode unfiltered`. This route validates the native opaque 16-bit PNG and records that no HDR radiometric bridge, physical scene lighting or OIDN was evaluated. It is appropriate for construction/anatomy/constraint views, not material realism claims.

The default diagnostic colour mode is single gray. Set `workbench_color_type` to `MATERIAL` or `OBJECT` explicitly when authored fibre/state colours carry meaning. The choice is recorded in the job and adaptations. All source scene bytes remain unchanged.

For legible construction labels, hide only the declared source text objects with the existing `visible_prefixes` capture option, then add flat high-contrast labels in empty image bands. Keep the overlay mask and verify native pixels outside it remain exact.

## Guards and acceptance

All routes require frozen source/dependency hashes, a captured pipeline snapshot, fresh outputs, and cache validation of actual output bytes. No raw or denoised result automatically passes visual/noise/detail acceptance. Preserve rejected attempts and selected baselines. A change in the denoising guides must be checked against the retained raw beauty before attributing any image change to the material.


## Explicit two-component filtering

`--mode component-recombined --component-filter filter_both` filters the positive transmission/emission component and the physical diffuse/glossy/environment/volume remainder separately. Both use the declared bounded auxiliary features. The original tiny negative transmission pass terms and float64 reconstruction residual bypass filtering and are restored unchanged. Values are never discarded or clipped. The existing `2^-16` relative arithmetic guard remains fail-closed; it is an integration threshold, not material calibration.

The default policy remains `retain_transmission`, preserving the qualified LCD/CRT route and its exact arithmetic order. No failure switches policies automatically. A policy on another capture mode is rejected before output creation. Job, source, pipeline, both filter scopes and all retained linear outputs are verified by the cache guard; a filter-policy mismatch invalidates reuse.

Both policies retain original pass EXRs and identity/native16 checks. `filter_both` also retains the filtered transmission linear EXR and its review PNG. Exact float32 EXR serialization uses the already-tested `save_render` route; the prior `image.save` defect is not reintroduced. Pixel noise and material detail require explicit native review after these numeric gates. No capture automatically updates a preferred visual reference.
