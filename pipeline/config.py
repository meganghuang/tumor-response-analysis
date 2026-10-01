"""Configuration loading and path resolution."""
from __future__ import annotations

import pathlib
from typing import Any

import yaml

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = REPO_ROOT / "config" / "config.yaml"


class Config(dict):
    """Dict with attribute access and repo-relative path resolution."""

    def __getattr__(self, key: str) -> Any:
        try:
            value = self[key]
        except KeyError as exc:  # pragma: no cover - defensive
            raise AttributeError(key) from exc
        return Config(value) if isinstance(value, dict) else value

    def path(self, kind: str, *parts: str) -> pathlib.Path:
        """Resolve a configured directory (``raw``/``processed``/``results``)."""
        base = REPO_ROOT / self["paths"][kind]
        base.mkdir(parents=True, exist_ok=True)
        return base.joinpath(*parts) if parts else base


def load_config(path: str | pathlib.Path | None = None) -> Config:
    with open(path or DEFAULT_CONFIG, "r", encoding="utf-8") as fh:
        return Config(yaml.safe_load(fh))
