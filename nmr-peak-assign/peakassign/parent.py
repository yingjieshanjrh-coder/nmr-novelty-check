"""Access to the parent project (nmr-novelty-check).

This tool lives in a subfolder of nmr-novelty-check and borrows its shift
predictor and its SQLite index. Nothing in the parent folder is written to:
the package is imported as it is and the index is opened read-only.
"""
import os
import sys
from pathlib import Path

PARENT = Path(__file__).resolve().parents[2]


def ensure_importable():
    """Make `import nmrcheck` work even when the parent is not pip-installed."""
    try:
        import nmrcheck  # noqa: F401
    except ImportError:
        if not (PARENT / "nmrcheck" / "__init__.py").exists():
            raise ImportError(
                "the nmrcheck package was not found; keep nmr-peak-assign "
                f"inside the nmr-novelty-check folder (looked in {PARENT})")
        sys.path.insert(0, str(PARENT))


def default_data_dir():
    return Path(os.environ.get("NMRCHECK_DATA") or PARENT / "data")


def open_index(data_dir=None):
    """Return the parent's Index, or None when no index has been built."""
    ensure_importable()
    from nmrcheck.index import Index, index_path

    data_dir = Path(data_dir) if data_dir else default_data_dir()
    if not index_path(data_dir).exists():
        return None
    return Index(data_dir)
