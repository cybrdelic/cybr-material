"""Resolve the one namespaced generation owner without legacy module injection."""
from functools import lru_cache
import importlib
from pathlib import Path
import sys
from .core import ContractError


@lru_cache(maxsize=1)
def api():
    repository = Path(__file__).resolve().parents[4]
    local = repository / "cybr_material"
    existing = sys.modules.get("cybr_material")
    if local.is_dir() and existing is not None:
        if not Path(existing.__file__).resolve().is_relative_to(local.resolve()):
            raise ContractError("A different MATERIAL implementation is already imported")
    old_path = sys.path[:]
    try:
        if local.is_dir():
            sys.path.insert(0, str(repository))
        module = importlib.import_module("cybr_material.generation")
    except ModuleNotFoundError as error:
        if error.name == "cybr_material":
            raise ContractError("Install the cybr-material source package before using standalone recovery") from error
        raise
    finally:
        sys.path[:] = old_path
    if module.RECIPE_REVISION != "material-3.2-de4146ac":
        raise ContractError("Unsupported installed MATERIAL recipe revision")
    return module
