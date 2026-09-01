from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from prepare_snippy_runtime import validate_runtime


def test_runtime_must_be_named_child_of_tmp() -> None:
    assert validate_runtime(Path("/tmp/hpylori-snippy-runtime")) == Path("/tmp/hpylori-snippy-runtime")
    with pytest.raises(ValueError):
        validate_runtime(Path("/tmp"))
    with pytest.raises(ValueError):
        validate_runtime(Path("/tmp/has space"))
