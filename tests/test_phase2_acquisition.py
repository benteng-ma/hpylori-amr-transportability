from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from acquire_phase2 import ncbi_tasks, zenodo_tasks


def test_frozen_assembly_acquisition_set():
    tasks = ncbi_tasks(ROOT)
    assert len(tasks) == 474
    assert len({(task.dataset_id, task.isolate_id) for task in tasks}) == 474
    assert all(task.source_url.startswith("https://ftp.ncbi.nlm.nih.gov/") for task in tasks)


def test_frozen_zenodo_read_acquisition_set():
    raw_record = ROOT / "metadata/raw/zenodo_10369064_2026-08-30.json"
    if not raw_record.is_file():
        pytest.skip("Zenodo source metadata are intentionally excluded from the public release")
    tasks = zenodo_tasks(ROOT)
    assert len(tasks) == 52
    assert len({task.isolate_id for task in tasks}) == 52
    assert all(task.expected_bytes and task.expected_checksum for task in tasks)
