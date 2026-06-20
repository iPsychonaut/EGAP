#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
record_provenance.py

Capture and reload pipeline provenance: the external commands EGAP runs, the
stage each ran in, and the version of each program used.

Every external tool invocation flows through ``run_subprocess_cmd``
(:mod:`subprocess_runner`), which calls :func:`record_command` here. Each call
appends a structured entry to ``<output_dir>/<sample_id>_provenance.json``,
anchored off the active per-sample log file set by
``initialize_logging_environment``. The HTML reporter later reloads this with
:func:`load_provenance` to build the run-overview section.

For runs completed before provenance capture existed (no JSON on disk),
:func:`load_provenance` falls back to parsing the per-sample log's ``CMD:`` and
``-> Running <stage>:`` lines, so older runs still get a provenance section.

Program versions are resolved by :func:`resolve_versions`, which probes
``<tool> --version`` first and falls back to the conda environment metadata in
``$CONDA_PREFIX/conda-meta``.

Stage:
    Cross-cutting (provenance capture + report input)

Created on 2026-06-13

Author: Ian Bollinger (ian.bollinger@entheome.org / ian.michael.bollinger@gmail.com)
"""
import os
import re
import sys
import json
import glob
import subprocess
from datetime import datetime


# Generic shell/util tokens that are commands but not "programs" worth a version.
_SKIP_TOOLS = {
    "bash", "sh", "zsh", "env", "cd", "true", "false", "awk", "gawk", "sed",
    "grep", "egrep", "cat", "cp", "mv", "rm", "mkdir", "echo", "ln", "chmod",
    "tee", "sort", "uniq", "head", "tail", "cut", "tr", "xargs", "find",
}

# Map an invoked binary (lowercased, with .py/.sh stripped) to its conda package
# name when they differ. Anything not listed is looked up under its own name.
_CONDA_ALIASES = {
    "bbmerge": "bbmap", "reformat": "bbmap", "bbduk": "bbmap", "bbmap": "bbmap",
    "quast": "quast", "spades": "spades", "ragtag": "ragtag",
    "compleasm": "compleasm", "nanoplot": "nanoplot", "pilon": "pilon",
    "purge_dups": "purge_dups", "pbccs": "pbccs", "gfatools": "gfatools",
}


# --------------------------------------------------------------
# Capture
# --------------------------------------------------------------
def _provenance_path():
    """Return the provenance JSON path derived from the active log file.

    Reads ``log.DEFAULT_LOG_FILE`` (``<dir>/<sample_id>_log.txt``) at call time
    and returns the sibling ``<dir>/<sample_id>_provenance.json``. Returns
    ``None`` when logging has not been initialized.
    """
    try:
        import log
        log_file = getattr(log, "DEFAULT_LOG_FILE", None)
    except Exception:
        log_file = None
    if not log_file:
        return None
    if log_file.endswith("_log.txt"):
        return log_file[: -len("_log.txt")] + "_provenance.json"
    base, _ext = os.path.splitext(log_file)
    return base + "_provenance.json"


def _normalize_tool(token):
    """Reduce an invoked binary path to a bare tool name (basename, no .py/.sh)."""
    name = os.path.basename(str(token)).strip()
    for suffix in (".py", ".sh"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


def _derive_tool(cmd_list):
    """Best-effort program name from a command (list or shell string)."""
    if isinstance(cmd_list, str):
        parts = cmd_list.split()
        token = parts[0] if parts else ""
    elif cmd_list:
        token = cmd_list[0]
    else:
        token = ""
    return _normalize_tool(token)


def _infer_stage():
    """Infer the pipeline stage from the running script's name (sys.argv[0])."""
    argv0 = sys.argv[0] if sys.argv else ""
    stage = os.path.basename(argv0)
    if stage.endswith(".py"):
        stage = stage[:-3]
    return stage or "unknown"


def _append_entry(entry):
    """Append *entry* to the per-sample provenance JSON (no-op if unavailable)."""
    path = _provenance_path()
    if not path:
        return
    data = []
    if os.path.exists(path):
        try:
            with open(path, "r", encoding="utf-8") as fh:
                data = json.load(fh)
            if not isinstance(data, list):
                data = []
        except Exception:
            data = []
    data.append(entry)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=2)


def record_command(cmd_list, stage=None):
    """Record one executed command: log it (with version) and persist it.

    Writes a ``PROVENANCE:`` line to stdout (captured into the per-sample and
    session logs) so the command and its resolved program version survive even
    if the HTML report is never generated, and appends a structured entry to
    the per-sample provenance JSON. Never raises.

    Parameters
    ----------
    cmd_list : list of str or str
        The command as passed to ``run_subprocess_cmd``.
    stage : str, optional
        Pipeline stage; inferred from the running script name when omitted.
    """
    try:
        command = cmd_list if isinstance(cmd_list, str) else " ".join(str(c) for c in cmd_list)
        tool = _derive_tool(cmd_list)
        meaningful = _normalize_tool(tool).lower() not in _SKIP_TOOLS
        stage = stage or _infer_stage()
        version = _version_for(tool) if meaningful else ""
        # Step-by-step log line. run_subprocess_cmd already printed the full
        # command as a CMD line (which carries the input/output file paths);
        # this adds the resolved program version so both are durable in the log.
        if meaningful:
            vtxt = version if (version and version != "unknown") else "version unknown"
            print(f"PROVENANCE:\t[{stage}] {tool} ({vtxt})")
        _append_entry({
            "type": "command",
            "stage": stage,
            "tool": tool,
            "command": command,
            "version": version,
            "time": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        })
    except Exception:
        # Provenance must never break a pipeline step.
        return


def record_file(label, path, stage=None):
    """Record a produced/consumed artefact path: log it and persist it.

    Writes a ``PROVENANCE FILE:`` line to stdout (captured into the logs) so
    artefact locations are traceable even without the HTML report, and appends
    a file entry to the provenance JSON. Never raises.

    Parameters
    ----------
    label : str
        Human-readable role, e.g. ``"Flye assembly"`` or ``"Illumina forward (dedup)"``.
    path : str
        Filesystem path to the artefact.
    stage : str, optional
        Pipeline stage; inferred from the running script name when omitted.
    """
    try:
        exists = bool(path and os.path.exists(str(path)))
        disp = str(path) if path else "(none)"
        print(f"PROVENANCE FILE:\t{label}: {disp}" + ("" if exists else "  (not found yet)"))
        _append_entry({
            "type": "file",
            "stage": stage or _infer_stage(),
            "label": label,
            "path": os.path.abspath(str(path)) if path else "",
            "exists": exists,
            "time": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        })
    except Exception:
        return


# --------------------------------------------------------------
# Reload (manifest, with log-parse fallback)
# --------------------------------------------------------------
def _parse_log_for_commands(log_path):
    """Reconstruct ordered ``{stage, tool, command}`` entries from a per-sample log.

    Used when no provenance JSON exists (runs predating provenance capture).
    Tracks the current stage from ``-> Running <stage>:`` lines emitted by the
    orchestrator and collects subsequent ``CMD:\\t<command>`` lines.
    """
    if not log_path or not os.path.exists(log_path):
        return []
    entries = []
    current_stage = "unknown"
    run_re = re.compile(r"Running\s+([A-Za-z0-9_]+)\s*:")
    with open(log_path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            mr = run_re.search(line)
            if mr:
                current_stage = mr.group(1)
                continue
            stripped = line.strip()
            if stripped.startswith("CMD:"):
                command = stripped[len("CMD:"):].strip()
                if command:
                    entries.append({
                        "stage": current_stage,
                        "tool": _derive_tool(command),
                        "command": command,
                    })
    return entries


def load_provenance(output_dir, sample_id):
    """Load provenance entries for a sample.

    Prefers ``<output_dir>/<sample_id>_provenance.json``; falls back to parsing
    ``<output_dir>/<sample_id>_log.txt`` for runs without a manifest.

    Returns
    -------
    tuple of (list, str)
        ``(entries, source)`` where *entries* is a list of
        ``{stage, tool, command[, time]}`` dicts and *source* is ``"manifest"``,
        ``"log"``, or ``"none"``.
    """
    output_dir = os.path.abspath(output_dir)
    json_path = os.path.join(output_dir, f"{sample_id}_provenance.json")
    if os.path.exists(json_path):
        try:
            with open(json_path, "r", encoding="utf-8") as fh:
                data = json.load(fh)
            if isinstance(data, list) and data:
                return data, "manifest"
        except Exception:
            pass
    log_path = os.path.join(output_dir, f"{sample_id}_log.txt")
    entries = _parse_log_for_commands(log_path)
    if entries:
        return entries, "log"
    return [], "none"


def pipeline_steps(entries):
    """Return the ordered, de-duplicated list of stages from *entries*."""
    seen = []
    for e in entries:
        stage = e.get("stage") or "unknown"
        if stage not in seen:
            seen.append(stage)
    return seen


# --------------------------------------------------------------
# Versions: probe --version, fall back to conda-meta
# --------------------------------------------------------------
_FIRST_VERSIONISH = re.compile(r"\d+\.\d+")


def _first_version_line(text):
    """Return the first line containing a dotted number, trimmed."""
    for raw in (text or "").splitlines():
        line = raw.strip()
        if line and _FIRST_VERSIONISH.search(line):
            return line
    return None


def probe_version(tool):
    """Probe ``<tool> --version`` / ``-v`` and return the first version-like line.

    Returns ``None`` if the tool cannot be probed (not found, no version flag,
    timeout). Never raises.
    """
    import shutil
    if not shutil.which(tool):
        return None
    for flag in ("--version", "-v", "-V", "version"):
        try:
            r = subprocess.run([tool, flag], capture_output=True, text=True, timeout=10)
        except Exception:
            continue
        line = _first_version_line((r.stdout or "") + "\n" + (r.stderr or ""))
        if line:
            return line
    return None


_CONDA_CACHE = None


def conda_versions():
    """Return ``{package_name: version}`` parsed from ``$CONDA_PREFIX/conda-meta``.

    Cached after the first call. Returns ``{}`` when not in a conda environment.
    """
    global _CONDA_CACHE
    if _CONDA_CACHE is not None:
        return _CONDA_CACHE
    out = {}
    prefix = os.environ.get("CONDA_PREFIX")
    if prefix:
        for jf in glob.glob(os.path.join(prefix, "conda-meta", "*.json")):
            try:
                with open(jf, "r", encoding="utf-8") as fh:
                    meta = json.load(fh)
                name = meta.get("name")
                version = meta.get("version")
                if name and version:
                    out[name] = version
            except Exception:
                continue
    _CONDA_CACHE = out
    return out


def resolve_versions(tools):
    """Resolve a version string for each tool: probe first, then conda metadata.

    Parameters
    ----------
    tools : iterable of str
        Tool/binary names (as recorded in provenance).

    Returns
    -------
    dict
        ``{tool: version_string}``; value is ``"unknown"`` when neither source
        yields a version.
    """
    conda = conda_versions()
    versions = {}
    for tool in tools:
        norm = _normalize_tool(tool).lower()
        probed = probe_version(tool)
        if probed:
            versions[tool] = probed
            continue
        pkg = _CONDA_ALIASES.get(norm, norm)
        versions[tool] = conda.get(pkg) or conda.get(norm) or "unknown"
    return versions


_VERSION_CACHE = {}


def _version_for(tool):
    """Return a cached resolved version for *tool* (probe then conda)."""
    if tool in _VERSION_CACHE:
        return _VERSION_CACHE[tool]
    try:
        version = resolve_versions([tool]).get(tool, "unknown")
    except Exception:
        version = "unknown"
    _VERSION_CACHE[tool] = version
    return version


def versions_from_entries(entries):
    """Return ``{tool: version}`` captured at run time from command entries.

    Only includes entries with a known (non-empty, non-``"unknown"``) version,
    so the reporter can prefer run-time-captured versions over re-probing.
    """
    out = {}
    for e in entries:
        if e.get("type") == "file":
            continue
        tool = e.get("tool")
        version = e.get("version")
        if tool and version and version != "unknown" and tool not in out:
            out[tool] = version
    return out


def meaningful_tools(entries):
    """Return the ordered, de-duplicated meaningful program names from *entries*.

    Filters out generic shell/util tokens (see ``_SKIP_TOOLS``).
    """
    out = []
    for e in entries:
        if e.get("type") == "file":
            continue
        tool = e.get("tool") or ""
        if not tool or _normalize_tool(tool).lower() in _SKIP_TOOLS:
            continue
        if tool not in out:
            out.append(tool)
    return out
