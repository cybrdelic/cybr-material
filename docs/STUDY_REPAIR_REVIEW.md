# Fabrication and material study repair

This is an unreleased study branch based on v3.1.1, commit
`d8b35a52c80667f88c336c9150f64f82cefcf54f`. All seventeen original native
PNG renders are retained in `evidence/`; their display copies are in
`docs/images/v3_1_1/`. The new PNGs are real Blender Cycles renders, using
native 4K procedural maps and editable geometry. No photographic textures,
scans, image generation, sharpening, artificial grain or upscaling were used.

## What changed

| Outputs | Correction |
|---|---|
| Walnut and oak | Reduced contour-like atlas contrast and removed the diagrammatic knot. Added a meter-space 3D growth volume for Blender: independent saw planes per object, radial cell bundles, coherent end grain, earlywood-linked pores and oak ray variation. Vertical cabinet panels and turned timber use the appropriate fiber axis. |
| Marble | A 2.8 m slab with narrow primary seams, subsidiary fractures and restrained mineral intergrowth. Geometry carries thickness and eased edges. Kitchen top and waterfall mapping unfolds the same slab across the corner. |
| Travertine | Smaller irregular pores correlated with bedding, shallow cavity floors and no painted dark outlines. Honed vein-cut finish. |
| Brass and steel | Specific intact satin finishes, fine directional tooling and rare scores. The fabricated metal studies use freshly sheared edges; the steel is black oxide rather than random rust islands. |
| Porcelain | Clear glaze, fine firing flow, a shaped open bowl, actual wall thickness and an unglazed foot ring. Random sprinkled chips were removed. |
| Leather | Submillimeter irregular grain, relief below 55 microns, fine follicles and restrained dye variation. Flexible cut samples have thickness, curved edges and real thread geometry. Upholstery has closed padding, support structure and sewn perimeter seams. |
| Plaster | Broad overlapping trowel passes and leading lips replace high-frequency peppering and random delamination. |
| Linen | Wandering yarn coordinates, varied yarn widths, slubs, actual openings, physical thickness, a hem and localized resting folds. Upholstered weave has a backing. |
| Three still lifes | Shared corrected finishes, coherent timber face/end grain, shaped glaze, restrained metals and sewn hide/textile geometry. Original cameras and lighting are retained. |
| Walnut salon | A 5 mm woven rug with continuous bound edge and local corner lift replaces the corrugated sheet. Cushions have localized compression and seams. The chair has a seat cradle and back supports; stone edges are fabricated. Garden ground, planting and an open pergola establish the exterior. |
| Oak library | Separate pages, covers, rounded spines, raised binding bands, stamped title details and horizontal stacks. Desaturated binding colors vary. Uprights terminate at the shelf rather than projecting above it, with recessed wall brackets. |
| Stone kitchen | Real cabinet bodies, continuous fitted fronts, 4 mm reveals, recessed plinths, mounted pulls, 40 mm stone and a seating overhang above 400 mm. Added hood intake baffles, burner caps, pot supports and controls. |
| Parquet | Individual 3D saw cuts, finish variation, beveled joints and a routed 18 mm brass recess. The 3 mm brass strip is flush at 20 mm. Removed the raised comparison batten and extended the surrounding floor to eliminate presentation wedges. |

## Comparisons and evidence

`docs/images/comparisons/` contains all seventeen appearance comparisons.
The ten new material specimens deliberately show edges and curved surfaces,
so those pairs are **presentation comparisons**, not controlled shader tests.
The salon, library and still lifes retain their original cameras and light
positions. The kitchen's main camera is wider to include the refrigerator;
its original v3.1.1 camera remains in the editable scene for comparison.
Exterior geometry and physical material response have changed.

The interiors also include explicit use and construction details: a fitted
refrigerator/freezer with seals, hinges and ventilation, a formed sink with
drain and gooseneck faucet, nickel pulls and socket plates, and bronze trim.
Salon and library walls have distinct mineral paint finishes, reading lamps,
fitted storage or a writing console. Leather pads have front-face concavity,
seat compression and dark sewn welts following their rounded edges.

`path_traced/controlled_surfaces/` additionally uses the original ten close-up
cameras, lights and physical view widths. These matched surface proofs allow
anatomy and finish to be judged independently of the new specimen shapes.
Marble now samples a larger quarry slab, so its sampled mineral region changes.
JPEG contact sheets and comparisons are display derivatives. Native RGB16 PNGs
are the evidence; actual samples, devices, timings and hashes are recorded in
the manifests. OptiX was used locally on an RTX 4060; CPU remains supported.

## Reproduce locally

Use Blender 4.5.3 LTS and the Python dependencies in `source/requirements.txt`.
Run from the repository root; replace `OPTIX` with `CPU` when unavailable.

```sh
python source/generate_materials.py --resolution 4096
python source/pack_unity.py --pipeline hdrp
python source/pack_unity.py --pipeline urp
blender --factory-startup -b --python source/build_blender.py
blender --factory-startup -b --python source/render_path_traced.py -- --kind macros --device OPTIX --size 1024 --samples 768 --min-samples 128
blender --factory-startup -b --python source/render_path_traced.py -- --kind scenes --device OPTIX --size 1280 --samples 768 --min-samples 128
blender --factory-startup -b --python source/render_path_traced.py -- --kind architecture --device OPTIX --size 2048 --samples 512 --min-samples 128 --denoise
blender --factory-startup -b --python source/render_path_traced.py -- --kind macros --surface-proof --draft --use-final-maps --device OPTIX --size 1024 --samples 768 --min-samples 128
python source/finalize_evidence.py
python source/package_path_traced.py
python source/package_suite.py
python source/make_comparisons.py
python source/verify_maps.py
python source/verify_deliverable.py
```

Validate each of the three `blender/CYBR_Cycles_*.blend` files using
`blender --factory-startup -b <file> --python-exit-code 1 --python source/verify_path_traced.py`.
The architectural floor detail overrides the room settings to raw 768 samples.
Device availability is checked explicitly; a missing requested GPU fails the run.

## Honest limits

These are physical-scale procedural approximations, not measured reflectance
or anatomical scans. The Blender timber volume and exported finite wood atlases
are separate representations; the exported maps do not encode the 3D end-grain
volume or per-object saw plane. The wood remains an approximation and is not a
complete cellular simulation of a tree. The sparse gray/ochre marble is quieter
than many heavily veined Oro slabs; its presentation should not be described as
a measured match to a quarry sample.

Cloth and upholstery are deterministic authored resting shapes with localized
folds and compression, not cached cloth or pressure simulations. Small creases
from bending and contact are represented by geometry rather than a mesh-specific
wear bake. Open stone pores do not include undercut internal volumes. Books have
binding and page-block structure and abstract title marks, not authored readable
titles. Exterior planting is intentionally simple procedural geometry.

Validation checks native map dimensions and depth, normal conventions, exact
packed channels, intact metal coverage, physical scale, displacement pairing,
portable texture paths, book/shelf contacts, shelf termination, rug thickness,
kitchen overhang and stone thickness, parquet coverage and brass non-overlap.
These checks support correctness; they do not certify photographic realism.
