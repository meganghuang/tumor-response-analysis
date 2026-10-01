"""Download the raw expression matrix and metadata listed in config/config.yaml."""
from pathlib import Path

import requests

from tumor_response.config import load_config, project_path


def download(url: str, dest_dir: Path) -> Path:
    dest = dest_dir / url.rstrip("/").split("/")[-1]
    if dest.exists():
        print(f"Already downloaded: {dest.name}")
        return dest
    print(f"Downloading {url}")
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        with open(dest, "wb") as f:
            for chunk in r.iter_content(chunk_size=1 << 20):
                f.write(chunk)
    return dest


def main():
    cfg = load_config()["data"]
    raw_dir = project_path(cfg["raw_dir"])
    urls = [cfg["expression_url"], cfg["metadata_url"]]
    if not all(urls):
        raise SystemExit("Set data.expression_url and data.metadata_url in config/config.yaml first.")
    for url in urls:
        download(url, raw_dir)


if __name__ == "__main__":
    main()
