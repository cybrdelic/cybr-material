"""Small synthetic cache fixtures; no Blender process or production mutation."""
import copy
import hashlib
import json
import os
import pathlib
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
from cache_validation import PIPELINE_FILES, hash_file, validate_cache


class CacheValidationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        environment = mock.patch.dict(os.environ)
        environment.start(); self.addCleanup(environment.stop)
        os.environ.pop('MATERIAL_PATH_MAP_FILE', None)
        self.base = pathlib.Path(self.temp.name)
        self.root = self.base / 'cloud-eevee'; self.root.mkdir()
        for name in PIPELINE_FILES:
            (self.root / name).write_text('original pipeline ' + name)
        files = {name: hash_file(self.root / name) for name in PIPELINE_FILES}
        digest = hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest()
        self.snapshot = self.root / 'pipeline_snapshots' / digest
        self.snapshot.mkdir(parents=True)
        for name in files:
            (self.snapshot / name).write_bytes((self.root / name).read_bytes())
        source = self.root / 'source.blend'; source.write_bytes(b'original scene')
        self.texture = self.root / 'native-map.png'; self.texture.write_bytes(b'original map')
        raw = self.root / 'panel.png'; self.clean = self.root / 'CLEAN_OIDN_panel.png'
        guide_dir = self.root / 'panel_guides'; guide_dir.mkdir()
        self.guides = [guide_dir / (name + '_0007.exr') for name in ('beauty', 'albedo', 'normal')]
        outputs = [raw, self.clean, *self.guides, guide_dir / 'denoised_linear.exr']
        for path in outputs:
            path.write_bytes(('unchanged ' + path.name).encode())
        self.job = {'source': str(source), 'source_sha256': hash_file(source),
                    'output': str(raw), 'color_depth': 16}
        self.receipt = {
            'job': copy.deepcopy(self.job), 'source_sha256': hash_file(source),
            'pipeline': {'files': files, 'sha256': digest, 'snapshot': str(self.snapshot)},
            'image_dependencies': [{'path': str(self.texture), 'packed': False,
                                    'sha256': hash_file(self.texture)}],
            'raw_path': str(raw), 'clean_path': str(self.clean),
            'guide_paths': list(map(str, self.guides)), 'raw_guides_retained': True,
            'files_sha256': {str(path): hash_file(path) for path in outputs},
            'oidn_guide_contract': {'valid': True, 'violations': []},
            'color_bridge_native_precision': {'bit_depth': 16, 'opaque': True,
                                               'max_error_native': 0},
        }
        self.receipt_path = self.root / 'clean-receipt.json'

    def enable_relocation(self):
        helper = self.base / 'material-production/materials/path_config.py'
        helper.parent.mkdir(parents=True); helper.write_text('fixed relocation helper')
        config = self.base / 'path-map.json'
        config.write_text('{"schema_version":1,"prefixes":{}}')
        files = dict(self.receipt['pipeline']['files'])
        files['path_config.py'] = hash_file(helper)
        digest = hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest()
        snapshot = self.root / 'pipeline_snapshots' / digest; snapshot.mkdir()
        for name in files:
            source = helper if name == 'path_config.py' else self.root / name
            (snapshot / name).write_bytes(source.read_bytes())
        self.receipt['pipeline'] = {'files': files, 'sha256': digest, 'snapshot': str(snapshot)}
        self.receipt['relocation'] = {'config_path': str(config), 'config_sha256': hash_file(config),
                                      'helper_sha256': hash_file(helper), 'mappings': {}, 'images': []}
        os.environ['MATERIAL_PATH_MAP_FILE'] = str(config)
        return config, helper

    def check(self, expected_valid, message=''):
        self.receipt_path.write_text(json.dumps(self.receipt))
        result = validate_cache(self.job, self.receipt_path, self.root)
        self.assertEqual(result['valid'], expected_valid, result)
        if message:
            self.assertTrue(any(message in item for item in result['reasons']), result)

    def test_unchanged_capture_passes(self):
        self.check(True)

    def test_external_map_mutation_rejected(self):
        self.texture.write_bytes(b'new map, same scene')
        self.check(False, 'External image changed')

    def test_missing_external_map_rejected(self):
        self.texture.unlink()
        self.check(False, 'unreadable cache')

    def test_packed_image_is_bound_by_scene_not_external_path(self):
        self.receipt['image_dependencies'][0]['packed'] = True
        self.texture.unlink()
        self.check(True)

    def test_source_mutation_rejected(self):
        pathlib.Path(self.job['source']).write_bytes(b'modified scene')
        self.check(False, 'Source scene changed')

    def test_current_pipeline_change_rejected(self):
        (self.root / 'oidn_image_bridge.py').write_text('new pipeline')
        self.check(False, 'Current pipeline changed')

    def test_snapshot_corruption_rejected(self):
        (self.snapshot / 'render.py').write_text('corrupt frozen snapshot')
        self.check(False, 'Modified pipeline snapshot')

    def test_clean_png_mutation_rejected(self):
        self.clean.write_bytes(b'wrong preview')
        self.check(False, 'Capture output changed')

    def test_guide_mutation_rejected(self):
        self.guides[1].write_bytes(b'wrong guide')
        self.check(False, 'Capture output changed')

    def test_output_hash_omission_rejected(self):
        del self.receipt['files_sha256'][str(self.clean)]
        self.check(False, 'Required output hash missing')

    def test_invalid_guides_rejected(self):
        self.receipt['oidn_guide_contract']['valid'] = False
        self.check(False, 'Invalid OIDN guide contract')

    def test_truncated_png_precision_rejected(self):
        self.receipt['color_bridge_native_precision']['bit_depth'] = 8
        self.check(False, 'Native PNG bit depth')

    def test_changed_job_rejected(self):
        self.job['samples'] = 256
        self.check(False, 'Job differs')

    def test_legacy_incomplete_receipt_rejected(self):
        del self.receipt['image_dependencies']
        self.check(False, 'Incomplete')

    def test_malformed_receipt_rejected(self):
        self.receipt_path.write_text('{partial')
        result = validate_cache(self.job, self.receipt_path, self.root)
        self.assertFalse(result['valid'])

    def test_unchanged_relocation_context_passes(self):
        self.enable_relocation(); self.check(True)

    def test_changed_relocation_config_rejected(self):
        config, _ = self.enable_relocation(); config.write_text('changed config')
        self.check(False, 'configuration content changed')

    def test_removed_relocation_context_rejected(self):
        self.enable_relocation(); os.environ.pop('MATERIAL_PATH_MAP_FILE')
        self.check(False, 'configuration context changed')

    def test_new_relocation_context_rejects_old_unmapped_capture(self):
        os.environ['MATERIAL_PATH_MAP_FILE'] = str(self.base / 'new-map.json')
        self.check(False, 'configuration context changed')

    def test_changed_relocation_helper_rejected(self):
        _, helper = self.enable_relocation(); helper.write_text('new helper')
        self.check(False, 'Current pipeline changed: path_config.py')

    def test_relocated_config_path_rejected_even_with_same_bytes(self):
        config, _ = self.enable_relocation(); other = self.base / 'other-map.json'
        other.write_bytes(config.read_bytes()); os.environ['MATERIAL_PATH_MAP_FILE'] = str(other)
        self.check(False, 'configuration location changed')


if __name__ == '__main__':
    unittest.main(verbosity=2)
