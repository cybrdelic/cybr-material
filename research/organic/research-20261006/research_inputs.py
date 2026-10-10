"""Resolve explicit external artifacts without historical machine defaults."""
from pathlib import Path
import sys
_CACHE = {}
def required_input(name):
    if name in _CACHE:
        return _CACHE[name]
    flag = "--" + name
    if flag not in sys.argv:
        raise RuntimeError(f"Missing required research input: {flag} PATH. See README.md; no archived state or scene is bundled.")
    index = sys.argv.index(flag)
    if index + 1 >= len(sys.argv) or sys.argv[index + 1].startswith("--"):
        raise RuntimeError(f"{flag} requires a file path")
    path = Path(sys.argv[index + 1]).expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Research input does not exist: {path}")
    del sys.argv[index:index + 2]
    _CACHE[name] = path
    return path
