from pathlib import Path

import pytest

from run_checkm2 import resolve_database


def test_resolve_database_accepts_nested_download_directory(tmp_path: Path) -> None:
    database = tmp_path / "CheckM2_database" / "uniref100.KO.1.dmnd"
    database.parent.mkdir()
    database.write_bytes(b"test")
    assert resolve_database(tmp_path) == database.resolve()


def test_resolve_database_rejects_missing_database(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        resolve_database(tmp_path)
