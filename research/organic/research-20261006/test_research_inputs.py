"""Tests only for external-input argument handling; no research process is run."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import research_inputs

class InputTests(unittest.TestCase):
    def setUp(self):
        research_inputs._CACHE.clear()
    def test_missing_argument_fails(self):
        with patch('sys.argv', ['study.py']):
            with self.assertRaisesRegex(RuntimeError, 'baseline-scene'):
                research_inputs.required_input('baseline-scene')
    def test_missing_value_fails(self):
        with patch('sys.argv', ['study.py', '--baseline-scene']):
            with self.assertRaises(RuntimeError):
                research_inputs.required_input('baseline-scene')
    def test_missing_file_fails(self):
        with tempfile.TemporaryDirectory() as d:
            with patch('sys.argv', ['study.py', '--baseline-scene', str(Path(d)/'missing.blend')]):
                with self.assertRaises(FileNotFoundError):
                    research_inputs.required_input('baseline-scene')
    def test_file_path_and_own_argument_removal(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'scene.blend';p.write_bytes(b'fixture')
            args=['blender', '--background', '--', '--baseline-scene', str(p), '--other', '12']
            with patch('sys.argv', args):
                self.assertEqual(research_inputs.required_input('baseline-scene'),p.resolve())
                self.assertEqual(args,['blender','--background','--','--other','12'])
                self.assertEqual(research_inputs.required_input('baseline-scene'),p.resolve())
    def test_directory_is_not_artifact(self):
        with tempfile.TemporaryDirectory() as d:
            with patch('sys.argv',['study.py','--capture-profile',d]):
                with self.assertRaises(FileNotFoundError):
                    research_inputs.required_input('capture-profile')

    def test_public_contract_omits_recovery_inventory(self):
        root=Path(__file__).resolve().parent
        missing=json.loads((root/'MISSING_INPUTS.json').read_text())
        self.assertNotIn('entries',missing)
        self.assertEqual({x['flag'] for x in missing['external_file_inputs']},{'--baseline-scene','--capture-profile','--studio-source'})
        text=json.dumps(missing)
        for key in ('relative_name','sha256','bytes','shape'):
            self.assertNotIn('"'+key+'":',text)
        sources=json.loads((root/'SOURCE_MANIFEST.json').read_text())
        self.assertNotIn('input_archives',sources)
        self.assertTrue(all('original_sha256' not in row for row in sources['sources']))

if __name__=='__main__':
    unittest.main()
