"""Locate the repository package for legacy source-script entry points."""
from importlib import import_module
from pathlib import Path
import sys

_PACKAGE_ROOT = Path(__file__).resolve().parents[1]


def ensure_package_path():
    """Make this checkout's package importable without an installation step."""
    root = str(_PACKAGE_ROOT)
    if root not in sys.path:
        sys.path.insert(0, root)
    return _PACKAGE_ROOT


def expose_helper(namespace, name):
    """Alias an old module to its sole maintained implementation.

    Copy attribute references for callers using spec_from_file_location; normal
    imports receive the canonical module itself, including private attributes.
    """
    ensure_package_path()
    module = import_module("cybr_material." + name)
    expected = (_PACKAGE_ROOT / "cybr_material" / (name + ".py")).resolve()
    if Path(module.__file__).resolve() != expected:
        raise ImportError("MATERIAL helper resolved outside this checkout: " + name)
    legacy_name = namespace["__name__"]
    namespace.update({key: value for key, value in vars(module).items()
                      if not key.startswith("__")})
    sys.modules[legacy_name] = module
    return module
