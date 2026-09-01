#!/usr/bin/env python3
"""Download one public file with retries, provenance, and SHA-256."""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--retries", type=int, default=4)
    args = parser.parse_args()

    error: Exception | None = None
    for attempt in range(1, args.retries + 1):
        try:
            request = Request(
                args.url,
                headers={"User-Agent": "hpylori-amr-transportability/0.0.1"},
            )
            with urlopen(request, timeout=300) as response:
                args.output.parent.mkdir(parents=True, exist_ok=True)
                digest = hashlib.sha256()
                byte_count = 0
                with args.output.open("wb") as handle:
                    while chunk := response.read(1024 * 1024):
                        if not chunk:
                            continue
                        handle.write(chunk)
                        digest.update(chunk)
                        byte_count += len(chunk)
                metadata = {
                    "source_url": args.url,
                    "final_url": response.geturl(),
                    "retrieved_at": datetime.now(timezone.utc).isoformat(),
                    "http_status": response.status,
                    "content_type": response.headers.get_content_type(),
                    "etag": response.headers.get("etag"),
                    "last_modified": response.headers.get("last-modified"),
                    "bytes": byte_count,
                    "sha256": digest.hexdigest(),
                }
                args.output.with_suffix(args.output.suffix + ".meta.json").write_text(
                    json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
                )
                return
        except Exception as exc:  # retry only; final exception is raised below
            error = exc
            if attempt < args.retries:
                time.sleep(2 ** (attempt - 1))
    raise SystemExit(f"download failed after {args.retries} attempts: {error}")


if __name__ == "__main__":
    main()
