#!/usr/bin/env python3
"""Query NCBI Entrez and preserve raw machine-readable records."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen


BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
BATCH_SIZE = 200


def fetch(url: str) -> bytes:
    request = Request(url, headers={"User-Agent": "hpylori-amr-transportability/0.0.1"})
    with urlopen(request, timeout=180) as response:
        return response.read()


def batched(values: list[str], size: int = BATCH_SIZE):
    for start in range(0, len(values), size):
        yield values[start : start + size]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", required=True)
    parser.add_argument("--term", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--retmax", type=int, default=10000)
    args = parser.parse_args()

    search_url = f"{BASE}/esearch.fcgi?" + urlencode({
        "db": args.db, "term": args.term, "retmode": "json", "retmax": args.retmax
    })
    search_body = fetch(search_url)
    ids = json.loads(search_body)["esearchresult"]["idlist"]
    summary_urls: list[str] = []
    merged_summary: dict = {"result": {"uids": []}}
    for batch in batched(ids):
        summary_url = f"{BASE}/esummary.fcgi?" + urlencode({
            "db": args.db, "id": ",".join(batch), "retmode": "json"
        })
        summary_urls.append(summary_url)
        parsed = json.loads(fetch(summary_url))
        result = parsed.get("result", {})
        merged_summary["result"]["uids"].extend(result.get("uids", []))
        for uid in result.get("uids", []):
            if uid in result:
                merged_summary["result"][uid] = result[uid]

    payload = {
        "query": args.term,
        "database": args.db,
        "search_url": search_url,
        "summary_urls": summary_urls,
        "retrieved_at": datetime.now(timezone.utc).isoformat(),
        "ids": ids,
        "search": json.loads(search_body),
        "summary": merged_summary,
    }
    rendered = (json.dumps(payload, indent=2, ensure_ascii=False) + "\n").encode()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(rendered)
    args.output.with_suffix(args.output.suffix + ".sha256").write_text(
        hashlib.sha256(rendered).hexdigest() + "  " + args.output.name + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
