# Shared MATERIAL generation: implementation checkpoint

The new `cybr_material` package owns the live procedural helper equations,
versioned material recipes, metric normal construction and exchange encoding.
It does not own LIGHT's BSDF evaluation, sampling, PDFs or light transport.

Source baseline: `365bd17aed0e3bbe9d2e7896f2bfdce1b9b54a92`.
Recipe compatibility revision: `material-3.2-de4146ac`, corresponding to the
unchanged historical source lock at `de4146ac759164e29323b152d37f64417c00c3fb`.
The historical generator blob is `12406a4d2359b826c3ee3fcf412caabb0441e489`.

## Current adoption

The five old `source/` helper paths are compatibility shims into the package;
they contain no copied numerical or shader-node implementations. This also
preserves Blender's `timber_volume_shader.install` dependency closure: importing
that helper does not require NumPy or SciPy.

The new package entry point is available:

```powershell
python -m cybr_material.generation --output-root <fresh-directory>
```

Production map generation requires native 4096 resolution. Explicit diagnostic
evaluation can use smaller native sizes, without upscaling exchange maps:

```powershell
python -m cybr_material.generation --resolution 32 --draft --only 05_champagne_brass --output-root <fresh-directory>
```

**Caller migration remains blocked.** The existing root
`source/generate_materials.py` CLI and the recovery/native-parity consumers have
not been redirected to this API. The immutable vendor snapshot is therefore
still reached by existing recovery code. This checkpoint must not be described
as completed cross-consumer consolidation.

## Supported computation contract

`evaluate_fields(identifier, resolution, revision=..., diagnostic=...)` validates
identity, recipe revision and resolution before generation. It rejects nonfinite
raw fields before applying the original authored range clamps. Normalized output
channels are read-only and bound to implementation/channel hashes.

`encode_exchange_maps(fields)` revalidates the complete field contract and yields
the original RGB8/gray8/height16 integer channels. It preserves sample centers,
filters, finite-cut versus periodic derivatives, height/micro-height pairing,
metric scale, linen opacity and quantized DirectX green inversion.

`export_material` and `export_suite` require fresh destinations. Completion
receipts are written only after complete finite map encoding. Their status
describes computed procedural-map validation; `physical_qualification` and
`visual_acceptance` remain false. These authored recipes are not measured scans.

All historical source-lock values and expected selected-map hashes remain
unchanged. LF attributes preserve those byte identities on fresh Windows
checkouts; normalizing a lock or accepting a changed expected hash is not used.

## Observed checks

```powershell
python -m unittest discover -s tests -p test_material_generation.py -v
```

Windows/Python 3.12: **10 tests pass**. All ten original recipes at 32 pixels
produce all twelve PNG maps byte-identically to the unchanged historical CLI;
decoded integer arrays and material metadata are also identical. Tests reject
unknown IDs/revisions, invalid production resolution, output collisions,
nonfinite raw fields and changed channels. DirectX byte inversion is checked
without aliasing the OpenGL result.

The comparison runs the historical source unchanged in a temporary test root;
it does not modify the retained snapshot or original output maps. Small numeric
test maps are not scene-render deliverables.

## Outstanding acceptance

- Redirect the actual root CLI and all brass recovery/native-parity callers.
- Compare all ten native 4096 map sets and codecs, with measured resource use.
- Rerender all affected existing material macro, still-life and architecture
  examples at their requested settings, plus the agreed new example.
- Execute the new package checks on Linux and macOS; current results are Windows.
- Record final exact source refs and standardized before/after line accounting.

No full map set or scene render was generated for this checkpoint. No remote
branch, PR or release was published.
