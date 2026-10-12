# Opt-in granular matrix: DU5

This is a reusable original procedural material component. It does not change
any selected material, the ten-material suite, the 35-material audit registry,
or their protected source recipes. It does not add a new standalone repository.

`source/granular_matrix.py` owns the aggregate height, scalar mineral-pigment
response and approximate unresolved-roughness response. The one authoritative
parameter set is in `source/granular_matrix_profiles.json`. Consumers provide
world position in metres, pixel footprint in metres and a relief gain. The
helper can layer its outputs over an existing material without replacing its
chromatic endpoints or original normal chain.

SCENES owns the Victorville depositional fields, finite-stock accounting,
collar placement, plant-base geometry and CN palette binding. Those are not
hard-coded into this reusable material. Its preparation copies the pinned
MATERIAL code/profile into a generated runtime snapshot, never a second
independently authored recipe.

## Scope and evidence

DU5 produced a useful improvement over the retained DR matrix in an actual
native Cycles crop. It is not a finished photoreal surface. The shrub trunk and
large clast still show synthetic form. Earlier DU4 cellular relief was rejected
and is not a profile/default in this patch.

The isolated metre-scale slab was rendered at 800 × 600, 32 raw Cycles samples,
with an identity object matrix and with an equivalent nonuniform matrix
(0.25, 4, 2). The world vertices were exactly identical. Decoded linear RGB
matched exactly, with zero maximum error. This verifies transform invariance
in Blender 4.3.2 for that fixture; it is not a claim of measured soil optics.

Characteristic wavelengths are 16, 8, 4, 2 and 0.8 mm. They describe an irregular
aggregate spectrum, not individually counted mineral diameters. All heights
are optical shading relief; they do not create sediment mass or replace
explicit grain inventory.

The footprint filter and roughness transfer remain approximate. In particular,
the roughness coefficients are a heuristic proxy, not an exact EWA integral or
calibrated normal-variance/BRDF model. Do not label them physically validated.

Run headless profile contracts from this repository:

    python3 -m unittest discover -s tests -p test_granular_matrix.py

The scene adapter and raw renderer evidence are maintained in CYBR SCENES,
`scenes/victorville-road/source/du5`. The accepted DR source remains separate.

The standalone public fixture needs only Blender and this MATERIAL checkout:

    blender --background --factory-startup -t 2 --python-exit-code 1 --python tests/render_granular_matrix_fixture.py -- /path/to/fresh-fixture

Its neutral one-metre slab at 640 × 480 and 16 raw Cycles samples also passed
exact decoded RGB equality under the same equivalent transforms. Use the
project resource guard when running it on a shared machine.
