from pathlib import Path
import sys,json,time,numpy as np
R=Path(__file__).resolve().parents[1];root=R.parents[1];sys.path[:0]=[str(R/'source'),str(R/'base_recipe/vendor')]
from reference_construction import packed_grains as before
from fractured_construction import packed_grains as after
bands=[(.0025,.0075,.40,6500),(.0011,.0035,.16,18000),(.00035,.0011,.13,22000)]
a=before(512,.25,1918,bands,.13);b=after(512,.25,1918,bands,.13);native=json.loads((R/'native4096/receipt.json').read_text())
r={'tested_resolution':512,'coarse_primitive_hash_before':a['primitive_sha256'],'coarse_primitive_hash_after':b['primitive_sha256'],'exact_coarse_seeded_positions_radii':a['primitive_sha256']==b['primitive_sha256'],'same_primitives_as_native4096':a['primitive_sha256']==native['fracture_model']['coarse_primitive_sha256'],'exact_coarse_coverage_mask':bool(np.array_equal(a['mask'],b['mask'])),'grading_before':a['stats'],'grading_after':b['stats'],'height_modified':not np.array_equal(a['height_m'],b['height_m'])}
assert r['exact_coarse_seeded_positions_radii'] and r['same_primitives_as_native4096'] and r['exact_coarse_coverage_mask'] and r['height_modified'];(R/'preservation_receipt.json').write_text(json.dumps(r,indent=2));print(json.dumps(r))
