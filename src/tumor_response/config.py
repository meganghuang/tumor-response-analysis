from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def load_config(path: str | Path = PROJECT_ROOT / "config" / "config.yaml") -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def project_path(relative: str) -> Path:
    """Resolve a config path relative to the project root and make sure it exists."""
    path = PROJECT_ROOT / relative
    path.mkdir(parents=True, exist_ok=True)
    return path
