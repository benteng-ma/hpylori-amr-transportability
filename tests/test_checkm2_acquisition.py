from pathlib import Path

from acquire_checkm2_database import digest


def test_digest_known_content(tmp_path: Path) -> None:
    path = tmp_path / "payload"
    path.write_bytes(b"checkm2")
    assert digest(path, "md5") == "ea6fbf002e3433bae0798fbe601fd602"
