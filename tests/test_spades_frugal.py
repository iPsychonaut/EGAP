"""bin/assemble_spades.py: --frugal below 8 threads, overridable for the A/B."""
import pytest

import assemble_spades as sp


@pytest.mark.parametrize("threads,expect", [(4, True), (7, True), ("7", True), (8, False), (16, False), ("x", False)])
def test_frugal_below_eight_threads(monkeypatch, threads, expect):
    monkeypatch.delenv("EGAP_SPADES_FRUGAL", raising=False)
    assert sp.spades_frugal(threads) is expect


@pytest.mark.parametrize("forced,threads,expect", [("1", 16, True), ("0", 4, False)])
def test_env_override(monkeypatch, forced, threads, expect):
    monkeypatch.setenv("EGAP_SPADES_FRUGAL", forced)
    assert sp.spades_frugal(threads) is expect
