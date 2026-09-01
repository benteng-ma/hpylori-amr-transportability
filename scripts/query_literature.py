#!/usr/bin/env python3
"""Reproducible Europe PMC query with raw JSON preservation."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--query", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--page-size", type=int, default=100)
    args = parser.parse_args()

    params = urlencode({
        "query": args.query,
        "format": "json",
        "pageSize": args.page_size,
        "resultType": "core",
    })
    url = f"{API}?{params}"
    request = Request(url, headers={"User-Agent": "hpylori-amr-transportability/0.0.1"})
    with urlopen(request, timeout=120) as response:
        body = response.read()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(body)
    metadata = {
        "query": args.query,
        "url": url,
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
        "sha256": hashlib.sha256(body).hexdigest(),
        "bytes": len(body),
        "result_count": len(json.loads(body).get("resultList", {}).get("result", [])),
    }
    args.output.with_suffix(args.output.suffix + ".meta.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()

