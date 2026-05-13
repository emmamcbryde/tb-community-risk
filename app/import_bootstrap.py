from __future__ import annotations

from pathlib import Path
import sys


def _path_inside(child: Path, parent: Path) -> bool:
    try:
        child.resolve().relative_to(parent.resolve())
        return True
    except ValueError:
        return False


def ensure_repo_root(repo_root: Path | None = None) -> Path:
    """Ensure this repository wins local imports such as engine.dynamic.

    Streamlit can keep modules alive across page changes. If some other module
    named `engine` is already loaded, Python may not re-resolve the repository's
    `engine/` package even when the repo root is on sys.path.
    """

    if repo_root is None:
        repo_root = Path(__file__).resolve().parents[1]

    repo_root = repo_root.resolve()
    repo_root_str = str(repo_root)

    # Put repo root first and remove duplicate later entries.
    cleaned = []
    for item in sys.path:
        if not item:
            cleaned.append(item)
            continue
        try:
            if Path(item).resolve() == repo_root:
                continue
        except Exception:
            pass
        cleaned.append(item)

    sys.path[:] = [repo_root_str] + cleaned

    # If a non-repo `engine` module is already cached, remove it and children.
    engine_mod = sys.modules.get("engine")
    if engine_mod is not None:
        engine_file = getattr(engine_mod, "__file__", None)
        engine_paths = list(getattr(engine_mod, "__path__", []))

        in_repo = False

        if engine_file:
            in_repo = _path_inside(Path(engine_file), repo_root / "engine")

        if engine_paths:
            in_repo = any(
                Path(p).resolve() == (repo_root / "engine").resolve()
                for p in engine_paths
            )

        if not in_repo:
            for name in list(sys.modules):
                if name == "engine" or name.startswith("engine."):
                    sys.modules.pop(name, None)

    return repo_root