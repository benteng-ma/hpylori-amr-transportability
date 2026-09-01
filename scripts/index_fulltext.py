#!/usr/bin/env python3
"""Index JATS full-text XML, supplement links, and accession statements."""

from __future__ import annotations

import argparse
import csv
import re
import xml.etree.ElementTree as ET
from pathlib import Path


XLINK = "{http://www.w3.org/1999/xlink}href"
ACCESSION_RE = re.compile(r"\b(?:PRJ[ENDA][A-Z]?\d+|SR[APRX]\d+|ER[APRX]\d+|GCA_\d+\.\d+|GCF_\d+\.\d+)\b")


def text_of(element: ET.Element | None) -> str:
    if element is None:
        return ""
    return " ".join(" ".join(element.itertext()).split())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--article-output", required=True, type=Path)
    parser.add_argument("--supplement-output", required=True, type=Path)
    args = parser.parse_args()

    articles: list[dict[str, str]] = []
    supplements: list[dict[str, str]] = []
    for path in sorted(args.input_dir.glob("*.xml")):
        root = ET.parse(path).getroot()
        title = text_of(root.find(".//article-title"))
        ids = {node.attrib.get("pub-id-type", "unknown"): text_of(node) for node in root.findall(".//article-id")}
        body_text = text_of(root.find(".//body"))
        accessions = sorted(set(ACCESSION_RE.findall(body_text)))
        availability_sections = []
        for section in root.findall(".//sec"):
            section_title = text_of(section.find("title")).lower()
            if "data availability" in section_title or "data sharing" in section_title:
                availability_sections.append(text_of(section))
        articles.append({
            "source_file": path.name,
            "pmcid": ids.get("pmc", path.stem.removeprefix("PMC")),
            "pmid": ids.get("pmid", ""),
            "doi": ids.get("doi", ""),
            "title": title,
            "accessions": ";".join(accessions),
            "data_availability": " || ".join(availability_sections),
        })
        for node in root.findall(".//supplementary-material") + root.findall(".//supplement"):
            hrefs = []
            if node.attrib.get(XLINK):
                hrefs.append(node.attrib[XLINK])
            for descendant in node.iter():
                if descendant.attrib.get(XLINK):
                    hrefs.append(descendant.attrib[XLINK])
            supplements.append({
                "source_file": path.name,
                "pmcid": ids.get("pmc", path.stem.removeprefix("PMC")),
                "doi": ids.get("doi", ""),
                "supplement_id": node.attrib.get("id", ""),
                "label": text_of(node.find("label")),
                "caption": text_of(node.find("caption")),
                "href": ";".join(dict.fromkeys(hrefs)),
            })

    args.article_output.parent.mkdir(parents=True, exist_ok=True)
    with args.article_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=articles[0].keys())
        writer.writeheader()
        writer.writerows(articles)
    with args.supplement_output.open("w", newline="", encoding="utf-8") as handle:
        fields = ["source_file", "pmcid", "doi", "supplement_id", "label", "caption", "href"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(supplements)


if __name__ == "__main__":
    main()
