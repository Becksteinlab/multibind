from pathlib import Path

from .multibind import Multibind

_PROJECT_ROOT = Path(__file__).resolve().parents[1]

try:
    from importlib.metadata import version as _package_version

    __version__ = _package_version("multibind")
except Exception:  # pragma: no cover
    __version__ = "0.0.0"

# Poetry editable installs can leave metadata at the placeholder; use versioningit from a dev env.
if __version__ == "0.0.0":
    try:
        from versioningit import get_version

        __version__ = get_version(_PROJECT_ROOT)
    except Exception:  # pragma: no cover
        pass
