# Procedural material pipeline

The deliverable is the generator, the physical material metadata and the
repeatable rendering workflow, as well as the generated maps. A bitmap is an
output of the material model. It is not the starting point for inventing the
other channels.

There are no input photographs, scans or image-generated textures. Random
seeds vary finite structures and process conditions. They do not directly
produce independent per-pixel color or bump.

## Wood: sample anatomy, then apply finish

`source/wood_anatomy.py` samples a board from a radial tree-growth volume. The
trunk runs along Y; the board samples an XZ cross-section. Pith position, saw
depth, saw slope and small trunk curvature determine radial position on the
face. That radius crosses nonuniform annual growth intervals. An optional
oblique branch cylinder changes the surrounding growth coordinates and
produces a sound knot intersection without a painted dark center.

The growth field supplies latewood density, annual pigment variation and
earlywood locations. Oak's large vessels concentrate in earlywood. Walnut's
smaller vessels occur through more of the ring. Individual projected lumen
windows have anatomical radius, finite length, tapered elliptical intersections,
fractional pixel coverage and rounded cavity floors. Their depth depends on
their radius and partial filling. Placement also follows cell-bundle regions
instead of a uniform sprinkle.

Medullary rays are finite tapered plates. Oak has larger visible plates;
walnut's fine rays contribute less. Cell bundles are sampled in radial and
longitudinal coordinates at three physical scales. They turn with the growth
volume and continue between vessel openings. Broad annual pigment bands do not
turn into height grooves.

The finish model applies contact polishing, shallow abrasive marks and a few
dry checks. Those features modify color, roughness, relief and wear together.
Geometry and finish produce height independently of the color output.

Edit `source/wood_profiles.json` to change a wood recipe:

| Parameter | Meaning |
|---|---|
| `extent_m` | Physical width and height of the generated board cut |
| `pith_x_m`, `cut_depth_m`, `cut_slope` | The cut through the growth volume |
| `ring_width_m` | Typical annual growth interval; individual years vary |
| `ring_porous` | Whether the large vessels concentrate in earlywood |
| `branches` | Optional cylindrical branch intersections; position, tilt and radius in meters |
| `vessel_radius_um`, `vessel_length_mm` | Anatomical opening sizes |
| `vessel_count` | Candidate population before anatomical placement |
| `ray_*` | Population and dimensions of finite medullary plates |
| `height_scale_m` | Physical encoding scale for the height map |
| `base_srgb`, `late_srgb`, `pore_srgb`, `ray_srgb` | Intrinsic material colors |
| `roughness` | Finish baseline |
| `seed` | Reproducible variation of the recipe |

The metadata reads physical extent, roughness and height scale from this same
recipe so the shader and generated normals use consistent units. The woods are
finite cuts. A continuous 3D model does not automatically make a square crop
periodic. Use an appropriate cut on a board rather than claiming the edges are
seamless.

This is a simplified anatomy model, not a biological growth simulation. It
does not yet generate arbitrary end-grain cuts, complete branch growth history
or cell-level optical scattering. The branch transition is an approximation.

## Other surface structures

`source/generate_materials.py` contains the other recipes. Shared utilities in
`surface_math.py` and `physical_features.py` construct periodic process fields,
finite deposits, fracture growth and overlapping applications.

Deposits use variable irregular polygonal rims and varied shallow floors. The
earlier shared angular sine profile produced recognizable flower-shaped pits;
that profile was removed during visual review.

- Marble uses irregular mineral seam widths, fractured margins, adjacent
  breccia fragments and shallow calcite structures.
- Travertine places torn cavities in porous sediment beds. Varied cavity
  floors produce depth; mineral residues supply intrinsic color differences.
- Metals distinguish tool marks, abrasion, bare metal and dielectric oxide.
  Corrosion fronts carry granular deposits and pits rather than a single
  flat color mask.
- Porcelain cracks grow and stop at earlier cracks. Glaze ripples, pinholes
  and chips have different physical dimensions.
- Leather combines stretched crease contours, follicle pits, compression,
  burnishing and grouped finish scuffs.
- Plaster layers tapered trowel passes, leading lips, polished interiors,
  mineral aggregate, pinholes and delamination.
- Linen constructs warp and weft with continuous over-under bends, fibrils,
  slubs and opacity derived from yarn coverage.

## Height and normals

The height field encodes meters as `(H - 0.5) × height_scale_m`. Normal gradients
use the same physical width and height. OpenGL and DirectX outputs differ only
in their green channel.

The full field is divided into a smooth geometry band and a residual band.
`Height_Macro.png` moves geometry. `Normal_Micro_OpenGL.png` supplies the
remaining fine relief. Combining geometry displacement with the full normal
would double-count relief. For finite wood cuts, the Blender shader derives
the selected height band's slope for the actual UV mapping.

## Run and verify

```bash
python -m pip install -r source/requirements.txt
python source/build_suite.py
```

For a wood iteration, generate the native maps and review a smaller Cycles
render before spending time on all final frames:

```bash
python source/generate_materials.py --resolution 4096 \
  --only 03_american_walnut 04_fumed_oak
blender -b --python source/render_path_traced.py -- \
  --kind macros --draft --use-final-maps \
  --only 03_american_walnut 04_fumed_oak \
  --size 768 --samples 192 --min-samples 64 --threshold .01
```

Drafts stay outside the packaged deliverable. Final scenes use Cycles, physical
dimensions, true displaced inspection patches and native 16-bit output. PNG map
validation checks dimensions, normal convention, packed channels and opacity;
scene validation checks renderer settings, displacement pairing and portable
texture paths. The complete archive receives SHA-256 checksums and a ZIP CRC
check.

Pinned dependency versions and stored seeds make map generation reproducible.
Render timing and metadata can vary across machines. Automated checks establish
technical consistency; they cannot certify photographic realism.
