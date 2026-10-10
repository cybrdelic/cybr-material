# Organic and textile research source preservation

This additive snapshot preserves the later October 6 bark, wool-spinning and
carpet research that was not directly present in the earlier draft PR. It does
not replace the selected fleece/carpet r5 or bark r8 appearances and does not
start, resume or promote an experiment. No native scene, render, solver state,
cache, private storage identity, operator note or command receipt is included.

## Contents

- `bark/`: cell-wall and two-cell mechanics, constitutive and native-port work,
  growth surrogate, notch reference, periodic seam and conditional opening.
- `textile/`: short-yarn scaling, suspended spinning, contact migration, finite
  and generalized friction, flexible-friction bridge and process studies.
- `carpet/work/`: packed and migrating fibre cores, clearance/contact audits,
  rounded yarn unit, physical-length continuation and the later r11 tangent
  correction. The later source supersedes overlapping copies from the earlier
  carpet delta. Distinct historical revisions remain labeled by their folders.

## Status and limits

These studies are held research, not finished material models. Bark's bounded
reduced-model comparisons do not establish mature-bark formation or appearance.
Textile's bounded friction result does not qualify a production yarn. Carpet's
isolated span-36 tangent improvement does not accept a whole prefix: adjacent
span/station-plane tests failed and contact clearance remains unproven. Selected
appearances are unchanged. CPU continuation remains held.

## External inputs and validation

`SOURCE_MANIFEST.json` records source hashes and the limited input-path adapters.
`MISSING_INPUTS.json` documents public input interfaces only. Historical recovery
filenames, artifact hashes, exact array inventories and archive metadata are
excluded. The consuming source defines each study's expected numerical state.
External state must be supplied separately; no state is bundled. Several scripts
also need sibling source, Blender, compiler tooling and a compatible data layout.
This snapshot makes no complete-replay or clean-clone runtime claim.

Historical machine-specific scene/profile/studio references were replaced by
explicit command-line inputs. Supply `--baseline-scene PATH`,
`--capture-profile PATH` or `--studio-source PATH` when the relevant script asks
for one. Relative paths resolve from the invoking working directory. The helper
removes only its own arguments before existing parsers run; missing files fail
clearly. No placeholder silently substitutes a scene or recipe.

Only syntax, JSON parsing, archive integrity, source inventory and path-input
helper tests were checked for this preservation step. No numerical experiment,
renderer, native solver or appearance test was executed.
