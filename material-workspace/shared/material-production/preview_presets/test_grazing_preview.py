import json
import math
from pathlib import Path
import tempfile
import unittest
from grazing_preview import (ROOT, assert_source_unchanged, defaults, load_preset,
                             plan_preset, source_guard)


class PresetTests(unittest.TestCase):
    def test_default_exact_reference(self):
        p, plan = load_preset(), plan_preset()
        self.assertEqual(defaults()['default_preview_preset'], 'cybr_grazing_warm_v1')
        self.assertEqual(plan['composition']['mode'], 'reference')
        self.assertEqual(plan['camera']['location_m'], p['camera']['location_m'])
        self.assertEqual(plan['camera']['rotation_euler'], p['camera']['rotation_euler'])
        self.assertEqual(plan['camera']['ortho_scale'], p['camera']['ortho_scale'])
        self.assertEqual(plan['key']['energy'], p['key']['energy'])
        self.assertFalse(defaults()['automatic_batch_rerender'])
        self.assertEqual(plan['cycles']['samples'], 128)
        self.assertEqual(plan['cycles']['seed'], 1810)
        self.assertFalse(plan['cycles']['use_adaptive_sampling'])
        self.assertFalse(plan['cycles']['use_light_tree'])
        self.assertFalse(plan['cycles']['use_denoising'])
        self.assertEqual(p['pipeline']['oidn_version'], '2.5.1')
        self.assertEqual(p['job_defaults']['capture_mode'], 'guided')
        self.assertEqual(p['render']['image_color_depth'], '16')
        self.assertEqual(p['render']['dither_intensity'], 0)

    def test_uniform_scaling_and_translation(self):
        base = plan_preset()
        pivot = (1, 2, 3)
        plan = plan_preset(specimen_width_m=.6, specimen_pivot_m=pivot)
        factor = 2.4
        for block in ('camera', 'key'):
            for i in range(3):
                self.assertAlmostEqual(plan[block]['location_m'][i] - pivot[i],
                                       (base[block]['location_m'][i] - [0, 0, .008][i]) * factor)
        self.assertAlmostEqual(plan['key']['size'], base['key']['size'] * factor)
        self.assertAlmostEqual(plan['key']['energy'], base['key']['energy'] * factor ** 2)
        self.assertEqual(plan['camera']['rotation_euler'], base['camera']['rotation_euler'])
        self.assertEqual(plan['key']['rotation_euler'], base['key']['rotation_euler'])
        self.assertEqual(plan['world'], base['world'])
        self.assertFalse(plan['source_geometry_rescaled'])

    def test_reference_camera_is_not_reaimed_at_pivot(self):
        plan = plan_preset()
        backward = [plan['camera']['matrix_world'][i][2] for i in range(3)]
        delta = [plan['camera']['location_m'][i] - plan['specimen_pivot_m'][i] for i in range(3)]
        distance = math.sqrt(sum(v*v for v in delta))
        # Reference camera points elsewhere; silently aiming it at the light target is wrong.
        self.assertGreater(sum(abs(backward[i] - delta[i]/distance) for i in range(3)), .05)

    def test_detail_changes_framing_only(self):
        ref = plan_preset()
        detail = plan_preset(composition='detail', crop_width_m=.052,
                             crop_center_m=(.01, .02, .008))
        self.assertEqual(ref['key'], detail['key'])
        self.assertEqual(ref['world'], detail['world'])
        self.assertEqual(ref['camera']['rotation_euler'], detail['camera']['rotation_euler'])
        self.assertEqual(detail['camera']['ortho_scale'], .052)

    def test_full_fits_every_bound_corner(self):
        import itertools
        for width, height in [(640,640), (960,640), (640,960)]:
            plan = plan_preset(composition='full', width=width, height=height,
                               specimen_bounds_m=[[-.125,-.125,0], [.125,.125,.04]])
            camera = plan['camera']
            horizontal = camera['ortho_scale'] * min(1, width/height)
            vertical = camera['ortho_scale'] * min(1, height/width)
            for point in itertools.product(*zip([-.125,-.125,0], [.125,.125,.04])):
                delta = [point[i] - camera['location_m'][i] for i in range(3)]
                for axis, span in ((0,horizontal), (1,vertical)):
                    projection = sum(delta[i] * camera['matrix_world'][i][axis] for i in range(3))
                    self.assertLess(abs(projection), span/2)

    def test_bad_input_fails(self):
        for options in ({'specimen_width_m': 0}, {'specimen_width_m': float('nan')},
                        {'composition':'unknown'}, {'composition':'reference','crop_width_m':.1},
                        {'composition':'full'}, {'specimen_pivot_m':(0,0,float('inf'))},
                        {'width':-1}, {'height':1.5}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                plan_preset(**options)

    def test_output_freshness_and_source_immutability(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root/'source.blend'
            source.write_bytes(b'immutable source')
            source_path, output, digest = source_guard(source, root/'new')
            self.assertFalse(output.exists())
            assert_source_unchanged(source_path, digest)
            with self.assertRaises(FileExistsError):
                source_guard(source, root)
            with self.assertRaises(FileExistsError):
                source_guard(source, source)
            link = root/'alias'
            link.symlink_to(source)
            with self.assertRaises(FileExistsError):
                source_guard(source, link)
            source.write_bytes(b'changed')
            with self.assertRaises(RuntimeError):
                assert_source_unchanged(source_path, digest)

    def test_no_render_in_entrypoint(self):
        text = (ROOT/'prepare_preview.py').read_text()
        self.assertNotIn('bpy.ops.render', text)
        self.assertNotIn('save_userpref', text)
        self.assertIn('--disable-autoexec', text)


if __name__ == '__main__':
    unittest.main()
