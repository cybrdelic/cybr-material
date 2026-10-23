#!/usr/bin/env python3
"""Compatibility CLI for the single maintained MATERIAL generation API.

Recipe and encoding implementations live in cybr_material. Historical source
and output hashes remain in the immutable recovery snapshot, not on this path.
Existing output folders are protected; use --output-root for a fresh export.
"""
from _material_imports import ensure_package_path

ROOT = ensure_package_path()

from cybr_material.generation import main, normal_from_height, seam
from cybr_material.recipes import GENERATORS, SPECS


if __name__ == "__main__":
    main()
