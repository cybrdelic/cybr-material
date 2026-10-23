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


class GenerationContracts(unittest.TestCase):
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
