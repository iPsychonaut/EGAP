import json

import pytest

from audit.prototypes import checkpoint as cp


def _outputs(tmp_path, n=2, size=100):
    outs = []
    for i in range(n):
        p = tmp_path / f"out_{i}.fasta"
        p.write_bytes(b"A" * size)
        outs.append(p)
    return outs


def test_marker_roundtrip_allows_skip(tmp_path):
    outs = _outputs(tmp_path)
    m = cp.write_done_marker(tmp_path / "step.done.json", outs, step="decontaminate_assembly")
    ok, reasons = cp.step_can_skip(m, outs)
    assert ok and reasons == []


def test_marker_refuses_to_write_for_missing_output(tmp_path):
    with pytest.raises(FileNotFoundError):
        cp.write_done_marker(tmp_path / "x.json", [tmp_path / "nope.fasta"])
    assert not (tmp_path / "x.json").exists()


def test_marker_refuses_to_write_for_empty_output(tmp_path):
    p = tmp_path / "empty.fasta"
    p.write_bytes(b"")
    with pytest.raises(ValueError):
        cp.write_done_marker(tmp_path / "x.json", [p])


def test_truncated_output_after_marker_forces_rerun(tmp_path):
    """The exact --resume poisoning case: step completed, output later
    truncated (disk full, killed compression). Old EGAP would SKIP."""
    outs = _outputs(tmp_path)
    m = cp.write_done_marker(tmp_path / "s.json", outs)
    outs[0].write_bytes(b"A" * 10)
    ok, reasons = cp.step_can_skip(m, outs)
    assert not ok and any("size changed" in r for r in reasons)


def test_deleted_output_forces_rerun(tmp_path):
    outs = _outputs(tmp_path)
    m = cp.write_done_marker(tmp_path / "s.json", outs)
    outs[1].unlink()
    ok, reasons = cp.step_can_skip(m, outs)
    assert not ok and any("missing" in r for r in reasons)


def test_hash_detects_same_size_corruption(tmp_path):
    outs = _outputs(tmp_path, n=1)
    m = cp.write_done_marker(tmp_path / "s.json", outs, hash_files=True)
    outs[0].write_bytes(b"C" * 100)          # same size, different bytes
    ok, reasons = cp.step_can_skip(m, outs)
    assert not ok and any("hash" in r for r in reasons)
    ok2, _ = cp.step_can_skip(m, outs, check_hash=False)
    assert ok2                                # opt-out for multi-GB files


def test_legacy_free_text_marker_is_not_trusted(tmp_path):
    m = tmp_path / "decontamination_done.txt"
    m.write_text("Input assembly   : /x\nSequences kept   : 7\n")
    ok, reasons = cp.step_can_skip(m)
    assert not ok and reasons == ["no valid marker"]


def test_marker_from_older_code_declaring_fewer_outputs(tmp_path):
    outs = _outputs(tmp_path, n=2)
    m = cp.write_done_marker(tmp_path / "s.json", outs[:1])
    ok, reasons = cp.step_can_skip(m, expected=outs)
    assert not ok and "does not declare" in reasons[0]


def test_marker_is_json_provenance(tmp_path):
    outs = _outputs(tmp_path, n=1)
    m = cp.write_done_marker(tmp_path / "s.json", outs, step="polish", extra={"racon_rounds": 2})
    doc = json.loads(m.read_text())
    assert doc["step"] == "polish" and doc["extra"]["racon_rounds"] == 2
    assert doc["outputs"][0]["size"] == 100
