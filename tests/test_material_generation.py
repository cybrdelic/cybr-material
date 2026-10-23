"""Bounded parity against unchanged historical CLI bytes, plus API contracts."""
from contextlib import contextmanager
from dataclasses import replace
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from types import MappingProxyType
import unittest
from unittest.mock import patch

import numpy as np
from PIL import Image

from cybr_material import generation as api

ROOT = Path(__file__).resolve().parents[1]
ORACLE = ROOT / "material-workspace/shared/material-production/vendor/original"
PRODUCTION = ORACLE.parents[1]


@contextmanager
def production_modules():
    """Load the actual recovery package without retaining a path injection."""
    prior = sys.path[:]
    try:
        sys.path.insert(0, str(PRODUCTION))
        from materials import generation, recover_native
        yield generation, recover_native
    finally:
        sys.path[:] = prior


class GenerationParity(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="material-parity-")
        cls.root = Path(cls.temp.name)
        oracle = cls.root / "oracle"
        (oracle / "source").mkdir(parents=True)
        lock = json.loads((ORACLE / "SOURCE_LOCK.json").read_text())
        for entry in lock["files"]:
            src = ORACLE / "source" / Path(entry["path"]).name
            if hashlib.sha256(src.read_bytes()).hexdigest() != entry["sha256"]:
                raise AssertionError("Historical oracle changed: " + src.name)
            shutil.copyfile(src, oracle / "source" / src.name)
        # The original CLI runs unchanged into its own temporary root. No
        # production assets or the retained historical snapshot are writable.
        result = subprocess.run([sys.executable, str(oracle / "source/generate_materials.py"),
                                 "--resolution", "32", "--draft"],
                                capture_output=True, text=True)
        if result.returncode:
            raise AssertionError(result.stdout + result.stderr)
        cls.old = oracle / "draft_materials"
        cls.new = cls.root / "shared"
        api.export_suite(cls.new, resolution=32, diagnostic=True)
        cls.cli = cls.root / "root-cli"
        result = subprocess.run([sys.executable, str(ROOT / "source/generate_materials.py"),
                                 "--resolution", "32", "--draft", "--only", "05_champagne_brass",
                                 "--output-root", str(cls.cli)], cwd=cls.root,
                                capture_output=True, text=True)
        if result.returncode:
            raise AssertionError(result.stdout + result.stderr)

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def test_all_ten_original_recipes_and_twelve_pngs_exact(self):
        self.assertEqual(len(api.recipes.SPECS), 10)
        self.assertEqual(len(api.MAP_TYPES), 12)
        for spec in api.recipes.SPECS:
            for name in api.MAP_TYPES:
                with self.subTest(material=spec["id"], channel=name):
                    old = self.old / spec["id"] / name
                    new = self.new / "draft_materials" / spec["id"] / name
                    self.assertEqual(old.read_bytes(), new.read_bytes())
                    with Image.open(old) as a, Image.open(new) as b:
                        np.testing.assert_array_equal(np.asarray(a), np.asarray(b))

    def test_metadata_keeps_metric_scales_and_channel_contracts(self):
        for spec in api.recipes.SPECS:
            old = json.loads((self.old / spec["id"] / "material.json").read_text())
            new = json.loads((self.new / "draft_materials" / spec["id"] / "material.json").read_text())
            self.assertEqual(old, new)

    def test_diagnostic_cannot_claim_qualification(self):
        report = json.loads((self.new / "draft_materials/generation.json").read_text())
        self.assertTrue(report["diagnostic"])
        self.assertFalse(report["physical_qualification"])
        self.assertFalse(report["visual_acceptance"])

    def test_root_cli_from_another_directory_uses_the_shared_encoder(self):
        folder = self.cli / "draft_materials/05_champagne_brass"
        for name in api.MAP_TYPES:
            self.assertEqual((folder / name).read_bytes(),
                             (self.old / "05_champagne_brass" / name).read_bytes())
        report = json.loads((folder / "generation.json").read_text())
        self.assertEqual(report["recipe_revision"], api.RECIPE_REVISION)
        self.assertEqual(report["implementation"], api.implementation_identity())
        self.assertTrue(report["diagnostic"])
        self.assertFalse(report["physical_qualification"])

    def test_brass_adapter_uses_shared_fields_without_executing_vendor(self):
        with production_modules() as (generation, _):
            with patch.object(generation, "load", side_effect=AssertionError("legacy engine ran")):
                fields, report = generation.evaluate("05_champagne_brass", 32)
        expected = api.evaluate_fields("05_champagne_brass", 32, diagnostic=True)
        for name, channel in expected.channels.items():
            np.testing.assert_array_equal(fields[name], channel)
        np.testing.assert_array_equal(fields["Normal_OpenGL"],
                                      api.normal_from_height(expected.channels["Height"], expected.spec, 32))
        self.assertEqual(report["recipe_revision"], api.RECIPE_REVISION)
        self.assertEqual(report["source_lock_sha256"], hashlib.sha256((ORACLE / "SOURCE_LOCK.json").read_bytes()).hexdigest())
        self.assertEqual(report["original_commit"], api.HISTORICAL_COMMIT)
        self.assertFalse(report["physical_qualification"])

    def test_native_recovery_dispatch_stops_before_full_resolution_allocation(self):
        """Wiring/preflight proof only: no native4096 recovery is claimed."""
        class StopBeforeAllocation(RuntimeError):
            pass
        with production_modules() as (_, recovery), tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "not-generated"
            before = sys.argv[:]
            with patch.object(api, "export_suite", side_effect=StopBeforeAllocation) as export:
                with self.assertRaises(StopBeforeAllocation):
                    recovery.recover_brass(output)
            export.assert_called_once_with(output, ["05_champagne_brass"], 4096,
                                           revision=api.RECIPE_REVISION)
            self.assertEqual(sys.argv, before)
            self.assertFalse(output.exists())

    def test_native_parity_caller_encodes_all_maps_in_a_bounded_wiring_harness(self):
        """Actual caller at a diagnostic substitution, not selected4096 proof."""
        with production_modules():
            from materials import native_parity
        test = self
        oracle = self.old / "05_champagne_brass"
        expected = api.evaluate_fields("05_champagne_brass", 32, diagnostic=True)

        class BoundedComparison:
            def __init__(self, identifier, out, compression):
                test.assertEqual(identifier, "05_champagne_brass")
                test.assertEqual(compression, 6)
                self.meta = dict(expected.spec)
                self.root, self.results = oracle, {}

            def compare(self, name, pixels, bits, encoded=False):
                test.assertTrue(encoded)
                test.assertEqual(bits, 16 if name in ("Height", "Height_Macro") else 8)
                with Image.open(self.root / (name + ".png")) as image:
                    np.testing.assert_array_equal(pixels, np.asarray(image))
                self.results[name] = True

            def snapshot(self, complete, extra):
                test.assertTrue(complete)
                test.assertEqual(len(self.results), 12)
                return extra

        with patch.object(native_parity, "Comparison", BoundedComparison), \
                patch.object(api, "evaluate_fields", return_value=expected) as evaluate:
            report = native_parity.brass("bounded-wiring-no-write")
        evaluate.assert_called_once_with("05_champagne_brass", 4096)
        self.assertEqual(report["recipe_revision"], api.RECIPE_REVISION)
        self.assertEqual(report["implementation"], api.implementation_identity())
        self.assertFalse(report["physical_qualification"])


class PackageInstallation(unittest.TestCase):
    def test_fresh_package_install_has_its_own_recipe_data_and_no_source_fallback(self):
        with tempfile.TemporaryDirectory(prefix="material-install-") as tmp:
            root = Path(tmp)
            build = root / "package"
            build.mkdir()
            shutil.copyfile(ROOT / "pyproject.toml", build / "pyproject.toml")
            shutil.copytree(ROOT / "cybr_material", build / "cybr_material",
                            ignore=shutil.ignore_patterns("__pycache__"))
            target = root / "installed"
            result = subprocess.run([sys.executable, "-m", "pip", "install", "--no-index", "--no-deps",
                                     "--no-build-isolation", "--disable-pip-version-check", "--target", str(target), str(build)],
                                    cwd=root, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            code = """import pathlib,sys
target=pathlib.Path(sys.argv[1]);sys.path.insert(0,str(target))
from cybr_material import generation as api
assert pathlib.Path(api.__file__).resolve().is_relative_to(target.resolve())
assert (pathlib.Path(api.__file__).parent/'wood_profiles.json').is_file()
fields=api.evaluate_fields('03_american_walnut',16,diagnostic=True)
assert fields.channels['BaseColor'].shape==(16,16,3)
assert len(list(api.encode_exchange_maps(fields)))==12
print(api.RECIPE_REVISION)
"""
            result = subprocess.run([sys.executable, "-I", "-c", code, str(target)],
                                    cwd=root, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(result.stdout.strip(), api.RECIPE_REVISION)


class GenerationContracts(unittest.TestCase):
    def test_row_hashing_preserves_c_order_digest_without_full_array_copy(self):
        value = np.arange(173 * 19 * 3, dtype=np.float32).reshape(173, 19, 3)
        for array in (value, value[::-1, ::2], value.transpose(1, 0, 2),
                      value.astype(np.uint16), np.asfortranarray(value)):
            self.assertEqual(api._array_sha256(array), hashlib.sha256(array.tobytes()).hexdigest())

    def test_diagnostic_intent_and_empty_requests_rejected_before_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "untouched"
            for flag in (1, "diagnostic", None):
                with self.assertRaises(api.ContractError):
                    api.export_suite(out, ["05_champagne_brass"], 32, diagnostic=flag)
                self.assertFalse(out.exists())
            for identifiers in ([], iter(())):
                with self.assertRaises(api.ContractError):
                    api.export_suite(out, identifiers, 32, diagnostic=True)
                self.assertFalse(out.exists())

    def test_source_changed_after_import_rejected_before_recipe_execution(self):
        changed = dict(api.implementation_identity(), **{"recipes.py": "changed-on-disk"})
        with patch.object(api, "implementation_identity", return_value=changed), \
                patch.object(api.recipes, "GENERATORS", None), self.assertRaises(api.ContractError):
            api.evaluate_fields("05_champagne_brass", 16, diagnostic=True)

    def test_unknown_identity_and_revision_rejected(self):
        for identifier, revision in (("missing", api.RECIPE_REVISION),
                                     ("05_champagne_brass", "unreviewed")):
            with self.assertRaises(api.ContractError):
                api.evaluate_fields(identifier, 32, revision=revision, diagnostic=True)

    def test_resolution_and_production_contract_rejected_before_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp) / "untouched"
            for n in (True, 0, 15, 4097, 32.0, float("nan"), 32):
                with self.subTest(n=n), self.assertRaises(api.ContractError):
                    api.export_suite(out, ["05_champagne_brass"], n)
                self.assertFalse(out.exists())
            with self.assertRaises(api.ContractError):
                api.export_suite(out, ["not-a-material"], 4096)
            self.assertFalse(out.exists())

    def test_output_collision_rejected_before_generation(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            (out / "draft_materials").mkdir()
            with patch.object(api, "evaluate_fields", side_effect=AssertionError("must not generate")):
                with self.assertRaises(FileExistsError):
                    api.export_suite(out, ["05_champagne_brass"], 32, diagnostic=True)

    def test_raw_nonfinite_is_rejected_before_clamping(self):
        n = 16
        bad = [np.zeros((n, n, 3), np.float32), *[np.zeros((n, n), np.float32) for _ in range(5)]]
        bad[1][0, 0] = np.nan
        generators = list(api.recipes.GENERATORS)
        generators[4] = lambda *_: bad
        with patch.object(api.recipes, "GENERATORS", generators), self.assertRaises(api.ContractError):
            api.evaluate_fields("05_champagne_brass", n, diagnostic=True)

    def test_authored_finite_out_of_range_fields_retain_legacy_clamp(self):
        n = 16
        raw = [np.full((n, n, 3), 2, np.float32), *[np.full((n, n), -1, np.float32) for _ in range(5)]]
        generators = list(api.recipes.GENERATORS)
        generators[4] = lambda *_: raw
        with patch.object(api.recipes, "GENERATORS", generators):
            fields = api.evaluate_fields("05_champagne_brass", n, diagnostic=True)
        self.assertEqual(float(fields.channels["BaseColor"].min()), 1)
        self.assertEqual(float(fields.channels["Height"].max()), 0)

    def test_mutated_or_nonfinite_channels_cannot_encode(self):
        fields = api.evaluate_fields("05_champagne_brass", 16, diagnostic=True)
        self.assertFalse(fields.channels["Height"].flags.writeable)
        for value in (np.full((16, 16), np.inf, np.float32), np.zeros((8, 8), np.float32),
                      np.zeros((16, 16), np.float32)):
            changed = dict(fields.channels)
            changed["Height"] = value
            with self.assertRaises(api.ContractError):
                list(api.encode_exchange_maps(replace(fields, channels=MappingProxyType(changed))))

    def test_directx_inverts_quantized_green_without_aliasing(self):
        fields = api.evaluate_fields("10_natural_linen", 33, diagnostic=True)
        maps = {name: pixels for name, pixels, _ in api.encode_exchange_maps(fields)}
        gl, dx = maps["Normal_OpenGL"], maps["Normal_DirectX"]
        np.testing.assert_array_equal(dx[..., 1], 255 - gl[..., 1])
        np.testing.assert_array_equal(dx[..., (0, 2)], gl[..., (0, 2)])
        self.assertFalse(np.shares_memory(gl, dx))


if __name__ == "__main__":
    unittest.main()
