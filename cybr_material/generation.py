"""One recipe/evaluation/encoding path for suite and recovery consumers.

Computed procedural maps are not measured-material or visual qualification.
Production maps retain the historical native 4096 resolution. Smaller runs
require an explicit diagnostic request and cannot be promoted to production.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import dataclass
import gc
import hashlib
import json
from pathlib import Path
import platform
from types import MappingProxyType
import zlib

import numpy as np
import PIL
from PIL import Image
import scipy
from scipy.ndimage import gaussian_filter

from . import recipes
from ._encoding_math import normal_from_height, seam
from .surface_math import clamp

RECIPE_REVISION = "material-3.2-de4146ac"
HISTORICAL_COMMIT = "de4146ac759164e29323b152d37f64417c00c3fb"
PRODUCTION_RESOLUTION = 4096
_CHANNELS = ("BaseColor", "Height", "Roughness", "AO", "Metallic", "WearMask", "Opacity")
MAP_TYPES = {
    "BaseColor.png": "sRGB RGB8; intrinsic color, no lighting",
    "Normal_OpenGL.png": "Non-Color tangent +Y RGB8",
    "Normal_DirectX.png": "Non-Color tangent -Y RGB8",
    "Normal_Micro_OpenGL.png": "Non-Color +Y RGB8; residual height only, paired with Height_Macro",
    "Height.png": "Non-Color linear grayscale16; displacement = (height - .5) * height_scale_m",
    "Height_Macro.png": "Non-Color linear grayscale16; band-limited geometry, same physical scale",
    "Roughness.png": "Non-Color linear grayscale8",
    "Metallic.png": "Non-Color grayscale8; exposed metal versus dielectric oxides",
    "AO.png": "Non-Color grayscale8; local cavity approximation, do not bake into color",
    "WearMask.png": "Non-Color grayscale8; material-specific abrasion or weathering",
    "ORM.png": "Non-Color RGB8; AO / roughness / metallic",
    "Opacity.png": "Non-Color grayscale8; woven openings for linen, opaque for other materials",
}


class ContractError(ValueError):
    """A request cannot produce maps under the supported recipe contract."""


def implementation_identity():
    root = Path(__file__).parent
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted([*root.glob("*.py"), root / "wood_profiles.json"])}


_IMPORTED_IMPLEMENTATION = implementation_identity()


def _checked_implementation():
    identity = implementation_identity()
    if identity != _IMPORTED_IMPLEMENTATION:
        raise ContractError("Generator source changed since import; start a fresh process")
    return identity


def runtime_identity():
    return {"python": platform.python_version(), "numpy": np.__version__,
            "scipy": scipy.__version__, "Pillow": PIL.__version__, "zlib": zlib.ZLIB_VERSION}


def _resolution(value, diagnostic):
    if type(diagnostic) is not bool:
        raise ContractError("Diagnostic intent must be an explicit boolean")
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise ContractError("Resolution must be an integer")
    n = int(value)
    if not 16 <= n <= PRODUCTION_RESOLUTION:
        raise ContractError("Supported native resolution is 16..4096")
    if not diagnostic and n != PRODUCTION_RESOLUTION:
        raise ContractError("Production exchange maps require native4096; request diagnostic explicitly")
    return n


def recipe(identifier, revision=RECIPE_REVISION):
    if revision != RECIPE_REVISION:
        raise ContractError("Unsupported recipe revision: " + str(revision))
    matches = [s for s in recipes.SPECS if s["id"] == identifier]
    if len(matches) != 1:
        raise ContractError("Unknown material identity: " + str(identifier))
    return deepcopy(matches[0])


@dataclass(frozen=True)
class MaterialFields:
    identifier: str
    revision: str
    resolution: int
    diagnostic: bool
    spec: object
    channels: object
    implementation: object
    channel_hashes: object


def _array_sha256(value):
    # Preserve the original C-order byte digest without a full native4K copy.
    # Sixty-four rows bound RGB float32 scratch to 3 MiB at resolution4096.
    digest = hashlib.sha256()
    for row in range(0, value.shape[0], 64):
        digest.update(value[row:row + 64].tobytes(order="C"))
    return digest.hexdigest()


def _hashes(channels):
    return {key: _array_sha256(value) for key, value in channels.items()}


def _validate_arrays(channels, n, *, normalized):
    if set(channels) != set(_CHANNELS):
        raise ContractError("The complete supported channel set is required")
    for key, value in channels.items():
        shape = (n, n, 3) if key == "BaseColor" else (n, n)
        if not isinstance(value, np.ndarray) or value.shape != shape or value.dtype.kind not in "fiu":
            raise ContractError("Unsupported channel shape/type: " + key)
        if not np.isfinite(value).all():
            raise ContractError("Nonfinite procedural channel: " + key)
        if normalized and (value.min() < 0 or value.max() > 1):
            raise ContractError("Out-of-range normalized channel: " + key)


def evaluate_fields(identifier, resolution=4096, *, revision=RECIPE_REVISION, diagnostic=False):
    n = _resolution(resolution, diagnostic)
    spec = recipe(identifier, revision)
    identity = _checked_implementation()
    index = next(i for i, value in enumerate(recipes.SPECS) if value["id"] == identifier)
    x = ((np.arange(n, dtype=np.float32) + .5) / n)[None, :]
    y = ((np.arange(n, dtype=np.float32) + .5) / n)[:, None]
    values = recipes.GENERATORS[index](x, y, n)
    if len(values) not in (6, 7):
        raise ContractError("Recipe returned an unsupported channel count")
    raw = dict(zip(_CHANNELS[:6], values[:6]))
    raw["Opacity"] = values[6] if len(values) == 7 else np.ones((n, n), np.float32)
    # Authored fields deliberately use the original range clamp. Invalid
    # arithmetic is rejected before that policy or any integer conversion.
    _validate_arrays(raw, n, normalized=False)
    channels = {key: clamp(raw[key]).astype(np.float32) for key in _CHANNELS[:6]}
    channels["Opacity"] = raw["Opacity"]
    mode = "nearest" if spec.get("tiling") == "finite-cut" else "wrap"
    channels["Height"] = gaussian_filter(channels["Height"], .35, mode=mode)
    _validate_arrays(channels, n, normalized=True)
    if implementation_identity() != identity:
        raise ContractError("Generator source changed during evaluation")
    hashes = _hashes(channels)
    for value in channels.values():
        value.setflags(write=False)
    return MaterialFields(identifier, revision, n, bool(diagnostic), MappingProxyType(spec),
                          MappingProxyType(channels), MappingProxyType(identity), MappingProxyType(hashes))


def validate_fields(fields):
    if not isinstance(fields, MaterialFields):
        raise ContractError("Encode fields returned by evaluate_fields")
    _resolution(fields.resolution, fields.diagnostic)
    if dict(fields.spec) != recipe(fields.identifier, fields.revision):
        raise ContractError("Recipe metadata changed")
    _validate_arrays(fields.channels, fields.resolution, normalized=True)
    if dict(fields.channel_hashes) != _hashes(fields.channels):
        raise ContractError("Generated channels changed after evaluation")
    if dict(fields.implementation) != _checked_implementation():
        raise ContractError("Generator source changed before encoding")


def encode_exchange_maps(fields):
    """Yield (channel name, native integer pixels, bit depth) once per map.

    A caller may stream these into PNG files or an exact native-pixel comparison.
    DirectX inversion operates on already-quantized OpenGL bytes, as before.
    """
    validate_fields(fields)
    c, s, n = fields.channels, fields.spec, fields.resolution
    h = c["Height"]
    mode = "nearest" if s.get("tiling") == "finite-cut" else "wrap"
    sigma = max(.5, n / 1024)
    if s.get("tiling") == "finite-cut":
        sigma = max(sigma, .62 * s["preview_diameter_m"] * 1.8 / 384 * n / s["tile_m"])
    macro = gaussian_filter(h, sigma, mode=mode)
    quantize = lambda a, bits=8: np.rint(clamp(a) * (65535 if bits == 16 else 255)).astype(np.uint16 if bits == 16 else np.uint8)
    yield "BaseColor", quantize(c["BaseColor"]), 8
    yield "Height", quantize(h, 16), 16
    yield "Height_Macro", quantize(macro, 16), 16
    yield "Normal_Micro_OpenGL", quantize(normal_from_height(h - macro, s, n)), 8
    for key in ("Roughness", "AO", "Metallic", "WearMask", "Opacity"):
        yield key, quantize(c[key]), 8
    yield "ORM", quantize(np.stack((c["AO"], c["Roughness"], c["Metallic"]), axis=-1)), 8
    normal = quantize(normal_from_height(h, s, n))
    yield "Normal_OpenGL", normal, 8
    directx = normal.copy()
    directx[..., 1] = 255 - directx[..., 1]
    yield "Normal_DirectX", directx, 8


def export_material(fields, destination):
    validate_fields(fields)
    out = Path(destination)
    out.mkdir(parents=True, exist_ok=False)
    encoded = {}
    for name, pixels, bits in encode_exchange_maps(fields):
        path = out / (name + ".png")
        Image.fromarray(pixels).save(path, compress_level=6)
        encoded[name] = {"bits": bits, "pixel_sha256": _array_sha256(pixels),
                         "png_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    if set(encoded) != {Path(name).stem for name in MAP_TYPES}:
        raise ContractError("Exchange-map coverage incomplete")
    spec = dict(fields.spec)
    (out / "material.json").write_text(json.dumps(dict(spec, resolution=fields.resolution, maps=MAP_TYPES), indent=2) + "\n", encoding="utf-8")
    size = min(768, fields.resolution)
    Image.fromarray(np.rint(fields.channels["BaseColor"] * 255).astype(np.uint8)).resize((size, size), Image.Resampling.LANCZOS).save(out / "swatch.jpg", quality=95)
    validate_fields(fields)
    report = {"status": "computed_exchange_maps_validated", "recipe_revision": fields.revision,
              "historical_commit": HISTORICAL_COMMIT, "diagnostic": fields.diagnostic,
              "resolution": [fields.resolution, fields.resolution], "units": "meters",
              "implementation": dict(fields.implementation), "runtime": runtime_identity(), "channels": encoded,
              "physical_qualification": False, "visual_acceptance": False,
              "scope": "Deterministic authored procedural maps; no measured-material or rendered-scene qualification"}
    with (out / "generation.json").open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2)
        stream.write("\n")
    return report


def export_suite(output_root, identifiers=None, resolution=4096, *, diagnostic=False, revision=RECIPE_REVISION):
    n = _resolution(resolution, diagnostic)
    ids = [s["id"] for s in recipes.SPECS] if identifiers is None else list(identifiers)
    if not ids:
        raise ContractError("At least one material identity is required")
    if len(ids) != len(set(ids)):
        raise ContractError("Duplicate material identities")
    specs = [recipe(identifier, revision) for identifier in ids]
    root = Path(output_root)
    folder = root / ("draft_materials" if diagnostic else "materials")
    folder.mkdir(parents=True, exist_ok=False)
    results, checks = {}, {}
    for spec in specs:
        print("GENERATING " + spec["id"] + " " + str(n) + "px", flush=True)
        fields = evaluate_fields(spec["id"], n, revision=revision, diagnostic=diagnostic)
        results[spec["id"]] = export_material(fields, folder / spec["id"])
        c = fields.channels
        checks[spec["id"]] = dict(height_min=float(c["Height"].min()), height_max=float(c["Height"].max()),
            relief_range_m=float((c["Height"].max() - c["Height"].min()) * spec["height_scale_m"]),
            roughness_min=float(c["Roughness"].min()), roughness_max=float(c["Roughness"].max()),
            metallic_min=float(c["Metallic"].min()), metallic_max=float(c["Metallic"].max()),
            base_color_seam=seam(c["BaseColor"]), height_seam=seam(c["Height"]), finite=True)
        del fields
        gc.collect()
    manifest = dict(name="CYBR MATERIAL / Fabrication and Physical Studies", version="3.2-study", resolution=n,
        workflow="metallic-roughness", normal="OpenGL +Y and DirectX -Y", base_color_space="sRGB",
        data_color_space="linear / Non-Color", height_bit_depth=16, height_midlevel=.5, units="meters",
        provenance="Ten original deterministic procedural materials. No source images, image generation, white-noise bump or upscaling. Not scans.",
        materials=deepcopy(recipes.SPECS))
    for filename, value in (("manifest.json", manifest), ("validation.json", checks), ("generation.json", {"status": "computed_suite_validated", "diagnostic": bool(diagnostic), "materials": results, "physical_qualification": False, "visual_acceptance": False})):
        with (folder / filename).open("x", encoding="utf-8") as stream:
            json.dump(value, stream, indent=2)
            stream.write("\n")
    print("MAPS_COMPLETE", flush=True)
    return results


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resolution", type=int, default=4096)
    parser.add_argument("--only", nargs="*", default=[])
    parser.add_argument("--draft", action="store_true", help="Explicit diagnostic maps; cannot promote to production")
    parser.add_argument("--output-root", type=Path, default=Path(__file__).resolve().parents[1], help="Fresh output root; existing map directories are protected")
    args = parser.parse_args(argv)
    export_suite(args.output_root, args.only or None, args.resolution, diagnostic=args.draft)


if __name__ == "__main__":
    main()
