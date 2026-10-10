# Rejected Victorville R2 · research provenance only

This visual candidate was rejected. Material 36 now selects the separately preserved R5D discrete-geometry ground; these historical R2 commands do not recreate that selection.

This second Victorville candidate responds to the request for much more dust, sand and dirt. It is a fully procedural, sediment-dominant regional interpretation. No scans, photographs or image-generation inputs are used. It is not an exact site reconstruction or a physically measured material.

## Physical change from the gravel-heavy first pass

The original 9,389 gravel/granule placements remain below the surface. A real deposited soil/sand mantle fills around and overtops them; the stones are anchored to the pre-deposition foundation instead of being lifted with the soil.

At native4096, the 60 × 60 cm repeat has:

- Exposed gravel: 2.9625% projected coverage, down from 33.4610%
- Dust/sand/dirt surface: 97.0375%
- Added sediment depth: 0.833–4.748 mm, mean 2.037 mm
- Maximum gravel height above the new soil: 3.813 mm
- 211 visibly exposed gravel bodies, with 792 low dirt agglomerates

The mantle varies with shallow drainage lows, broad sand deposits and sheltered/downwind zones near coarse stones. Dust coating increases near burial boundaries. Exposure, color, roughness and micro-relief share those same fields. Loose fine grains are concentrated in deposits rather than evenly decorating the surface. All percentages and dimensions are authored model outputs, not field measurements.

## User-selected preview lighting

The previews use the low-grazing asphalt disk-light rig requested by the user:

- Reference crop: 0.135 m
- Disk offset from surface focus: (−0.22, −0.16, +0.027) m
- Disk diameter: 0.055 m; power: 3.24 W; color: (1, 0.92, 0.82)
- For crop scale f, position and diameter scale by f; power scales by f²
- World color: (0.18, 0.21, 0.25), strength 0.04; no other lights
- Cycles CPU, two threads, 128 samples, seed 1810, adaptive sampling off
- Guided OpenImageDenoise 2.5.1; AgX / Medium High Contrast, exposure 0, gamma 1

Macro and detail views use 0.46 m and 0.13 m orthographic widths. The production images are 768 × 768; all image-based material channels are freshly evaluated at native4096. Low-angle lighting intentionally reveals millimetre and submillimetre relief; it is not a midday color-calibration view.

## Reproduce

Dependencies: Python with NumPy, SciPy and Pillow; Blender 4.3.2; an official CPU OpenImageDenoise installation for guided preview output. Set OIDN_BIN to its oidnDenoise executable. Large maps and Blender scenes are regenerated, not included in this compact source bundle.

    OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python source/generate_maps.py --resolution 4096 --out native4096
    OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python source/build_gravel_geometry.py native4096
    OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 python source/render_guided.py --view both --resolution 768 --samples 128 --oidn-bin "$OIDN_BIN"
    python source/test_ground.py

Adjust seed, tile_m, gravel_density, sediment_depth_m, sand_deposit_m, compaction_smoothing_m, loose_sand_scale, clod_density, dust_amount, crust_strength, color and height_scale_m. Physical feature sizes remain in metres. Height encoding is (Height − 0.5) × height_scale_m. The default height scale is 0.032 m.

The explicit gravel adapter clips seeded convex bodies from the same roof/edge planes, applies small upper-corner chips and 44–220 μm variable-width geometric edge bevels, and uses coupled mineral micro-relief. The soil receives its own metric normal and actual displaced mesh. The full-height normal is also exported for planar PBR workflows. Arbitrary deformed-mesh tangent transfer remains unqualified.

## Review status

This is a rejected visual candidate retained only for research provenance. The compact archive has no standalone native-output verifier; source/test_ground.py covers numerical contracts only. The first pass's large stones and coarse mineral markings were visibly stylized. Most are now buried/coated. The first dust-mantle version still appears soft/pillowy in places, with overly frequent dark fine-grain marks under grazing light. Further work should emphasize compacted/eroded flatter patches and less regular fine-scale contrast. This is not a photoreal or site-realism claim. No automatic selection or promotion is included.

## Regional provenance

- USDA NRCS Victorville series: predominantly granitic mixed alluvium, grayish-brown sandy loam and fine gravel. That named series concerns river terraces/floodplains and is context, not the assigned identity of this open-desert material: https://soilseries.sc.egov.usda.gov/OSD_Docs/V/VICTORVILLE.html
- USGS regional geology: sand, silt and gravel with differing felsic/granitic and local rock sources in the western Mojave/Victorville region: https://pubs.usgs.gov/publication/sir20235089/full
- NPS arid landforms: alluvial deposition, gravel pavements, wind and water transport provide the regional structural context: https://www.nps.gov/subjects/geology/arid-landforms.htm
- An openly accessible NPS Mojave road-edge photograph was visually inspected in the cloud browser. Pale gray-tan fines, mixed angular gravel and sorting were observed. It was not a Victorville site reference and was never used as an input texture: https://www.nps.gov/common/uploads/pwr/park/moja/673D14C1-F575-E2D5-CEFC0AF51F76CF60/673D14C1-F575-E2D5-CEFC0AF51F76CF60.jpg?autorotate=false&maxHeight=1020&maxWidth=1300
