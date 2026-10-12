# Opt-in longitudinal bark response

This original procedural closure layers fine height and scalar pigment variation
on an existing wood material. It preserves the source chromatic endpoints and
normal chain. No scan, photo texture or image generation is used. It changes no
selected MATERIAL recipe or catalogue entry.

The authoritative parameters live in `source/bark_flow_profiles.json`.
`source/bark_flow.py` reads that profile and consumes two point attributes:

- `dv_wood_flow_m`: continuous local coordinates across and along the branch,
  in prototype-local metres. Authored instance scaling also scales this spacing.
- `dv_wood_flow_valid`: one on the intended woody region and zero elsewhere.

Bump height coefficients are 0.70 and 0.19 millimetres in world-metre shader
units. These are field coefficients, not peak-to-peak measured bark heights.
The coarse stem form and junctions are actual geometry owned by SCENES; this
closure does not stand in for those shapes or choose the plant's placement.

This is a visual authored bark response, not measured plant optics. It does not
claim botanical calibration, complete scattering, or calibrated subpixel
normal-variance compensation. The particular scene still looks stylized.

Run the lightweight contracts from this repository:

    python3 -m unittest discover -s tests -p test_bark_flow.py

The coupled actual Cycles plant fixture is in CYBR SCENES under
`scenes/victorville-road/source/dv3/render_fixture.py`. It compares the retained
CM morphology with DV3 using the same camera and illumination. Refer to that
scene's evidence for render qualification; source tests alone do not establish
visual quality or soil/root collision correctness.
