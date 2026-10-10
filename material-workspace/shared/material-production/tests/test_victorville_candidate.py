"""Selected R5D geometry and rejected R2 research remain distinct."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
from materials.core import registry, material, validate_selection, ContractError, ROOT, sha
from materials.generation import evaluate, rejected_victorville_r2_fields

V = ROOT / 'vendor/victorville-r5d-20261009'
spec = importlib.util.spec_from_file_location('r5d_reconstruction_contract', V / 'reconstruct.py')
reconstruction = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reconstruction)

class VictorvilleCandidate(unittest.TestCase):
    def test_r5d_is_selected_with_exact_original_scene(self):
        self.assertEqual(len(registry()['materials']), 36)
        m = material('36')
        self.assertEqual(m['selected']['version'], 'approved_r5d_discrete_ground_20261009')
        review = json.loads((V / 'REVIEW.json').read_text())
        self.assertEqual(m['selected']['scene_sha256'], review['approved_scene_sha256'])
        self.assertEqual(m['pbr_contract']['representation'], 'discrete_geometry_plus_microbed')
        self.assertFalse(m['pbr_contract']['map_adapter_qualified'])
        self.assertFalse(m['pbr_contract']['seamless_tiling_qualified'])
        self.assertTrue(m['protected_appearance'])

    def test_selected_identity_never_falls_back_to_rejected_maps(self):
        with self.assertRaisesRegex(ContractError, 'discrete geometry'):
            evaluate('36_victorville_desert_ground', 64)

    def test_rejected_r2_remains_explicit_research(self):
        self.assertEqual(material('36')['rejected_research'][0]['status'], 'rejected_visual_r2')
        fields, report = rejected_victorville_r2_fields(32)
        self.assertEqual(fields['BaseColor'].shape, (32, 32, 3))
        self.assertIn('rejected visual R2', report['mode'])

    def test_original_approved_source_hashes(self):
        reconstruction.verify_source()
        self.assertEqual(sha(V/'src/build_diagnostic.py'), '320dd02e5ac661f4f9f3b470d614e0b3fb103614b0d2cd740aded5964ffa4246')
        self.assertEqual(sha(V/'src/fracture_prototypes.py'), '394bfc265be6e98cec9cd3f0bdbe8c2862a591bdfe0c0edd7262e3d92623c078')
        self.assertEqual(sha(V/'src/compacted_microbed.py'), '355576685bd6ef8492d0e4f2e5f7ab3963c58179a1335ec78ce924eecf378ec7')

    def test_original_geometry_receipt_is_explicit_history(self):
        receipt = json.loads((V/'receipts/geometry_d.json').read_text())
        self.assertEqual(receipt['particle_count'], 73733)
        self.assertEqual(receipt['extent_m'], .18)
        self.assertEqual(sum(row[1] for row in receipt['bands']), 73733)
        review = json.loads((V/'REVIEW.json').read_text())
        self.assertFalse(review['road_scene_approved'])
        self.assertFalse(review['scene_rebuild_byte_parity_qualified'])

    def test_reconstruction_refuses_overwrite_and_source_directory(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(reconstruction.subprocess, 'run') as command:
            with self.assertRaises(ValueError):
                reconstruction.reconstruct(directory, prepare_only=True)
            with self.assertRaises(ValueError):
                reconstruction.reconstruct(V/'never-write-here', prepare_only=True)
            command.assert_not_called()

    def test_reconstruct_inputs_in_fresh_directory_without_blender(self):
        with tempfile.TemporaryDirectory() as directory:
            result = reconstruction.reconstruct(Path(directory)/'fresh', blender='must-not-be-called', prepare_only=True)
            self.assertTrue(result['fracture_prototype_bytes_match'])
            self.assertTrue(result['microbed_array_contents_match'])
            self.assertFalse(result['geometry_rebuilt'])
            self.assertFalse(result['render_performed'])
            self.assertFalse(result['selected_scene_binary_replacement'])

if __name__ == '__main__':
    unittest.main()
