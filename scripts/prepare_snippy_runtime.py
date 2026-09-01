#!/usr/bin/env python3
"""Materialize Snippy helper scripts under a no-space Linux runtime path.

Snippy 4.6.0 interpolates its helper-script directory into shell commands
without quoting it.  The conda environment remains the software source of
truth; only the small Snippy script layer is copied to /tmp for execution.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_runtime(path: Path) -> Path:
    resolved = path.resolve()
    if resolved == Path("/tmp") or Path("/tmp") not in resolved.parents:
        raise ValueError("runtime must be a named directory below /tmp")
    if any(character.isspace() for character in str(resolved)):
        raise ValueError("runtime path must not contain whitespace")
    return resolved


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-prefix", type=Path, required=True)
    parser.add_argument("--runtime", type=Path, default=Path("/tmp/hpylori-snippy-runtime"))
    parser.add_argument("--manifest", type=Path, default=Path("results/logs/snippy_runtime_manifest.json"))
    args = parser.parse_args()

    source = args.source_prefix.resolve()
    runtime = validate_runtime(args.runtime)
    source_bin = source / "bin"
    scripts = sorted(source_bin.glob("snippy*"))
    snpeff_config = source / "etc/snpeff.config"
    if not scripts or not (source / "perl5").is_dir() or not snpeff_config.is_file():
        raise FileNotFoundError(f"invalid Snippy environment: {source}")

    if runtime.is_symlink():
        runtime.unlink()
    elif runtime.exists():
        shutil.rmtree(runtime)
    (runtime / "bin").mkdir(parents=True)
    (runtime / "etc").mkdir()
    copied = []
    for script in scripts:
        destination = runtime / "bin" / script.name
        shutil.copy2(script, destination)
        destination.chmod(destination.stat().st_mode | 0o111)
        copied.append({"name": script.name, "sha256": sha256(destination)})
    os.symlink(source / "perl5", runtime / "perl5", target_is_directory=True)
    shutil.copy2(snpeff_config, runtime / "etc/snpeff.config")

    path = f"{runtime / 'bin'}:{source_bin}:/usr/bin:/bin"
    version = subprocess.run(
        [str(runtime / "bin/snippy"), "--version"],
        env={**os.environ, "PATH": path}, capture_output=True, text=True, check=True,
    ).stdout.strip()
    manifest = args.manifest.resolve()
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps({
        "source_prefix": str(source), "runtime": str(runtime),
        "version": version, "copied_scripts": copied,
        "snpeff_config_sha256": sha256(runtime / "etc/snpeff.config"),
        "note": "Ephemeral no-space copy; source conda environment remains authoritative.",
    }, indent=2) + "\n", encoding="utf-8")
    print(version)


if __name__ == "__main__":
    main()
