"""EST_SIZE parsing and the report template's link names."""
import re
from pathlib import Path

import pytest

from estimate_runtime import parse_genome_size

REPO = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("text, expected", [
    ("43m", 43_000_000),
    ("1.2g", 1_200_000_000),
    ("5mb", 5_000_000),
    ("5 Mbp", 5_000_000),
    ("500k", 500_000),
    ("500kb", 500_000),
    ("25000000", None),   # no suffix: callers fall back to their default
    ("5x", None),         # unknown suffix: never guess a multiplier
    (None, None),
])
def test_parse_genome_size(text, expected):
    assert parse_genome_size(text) == expected


def test_template_links_are_supplied_by_reporter():
    """Every ``*_link`` name the summary template reads must be a render key."""
    template = (REPO / "resources" / "templates" / "EGAP_summary.html").read_text(encoding="utf-8")
    reporter = (REPO / "bin" / "html_reporter.py").read_text(encoding="utf-8")
    wanted = set(re.findall(r"\b(\w+_link)\b", template))
    missing = sorted(name for name in wanted if f'"{name}":' not in reporter)
    assert not missing, f"template reads link names the reporter never passes: {missing}"
