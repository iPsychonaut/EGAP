#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
sample_workbook.py

The sample input as a workbook (.xlsx): reading one, and generating one whose
every column header carries a hover tooltip.

Sheets
------
1. "Sample Information" (always the first sheet). The ground truth for which
   samples run. Same columns as the plain sample TSV.
2. "Sample Metadata". What the NCBI submission sheets need (BioSample,
   tissue, library strategy, instruments, iNaturalist observation).
3. "Sample Settings". RESERVED, not read yet; scope below.

Rule for every sheet after the first: a row's SAMPLE_ID must be a SAMPLE_ID of
sheet 1, and may appear once. A violation stops the run at start-up, because a
metadata or settings row that silently applies to no sample is worse than no
row. Sheets with no SAMPLE_ID column (notes, term lists) are ignored.

How the pipeline sees a workbook
--------------------------------
``workbook_to_tsv`` joins the metadata sheet onto the samples and writes
``<output_dir>/<workbook stem>_samples.tsv``. The orchestrator does this once
at start-up and hands that TSV to every stage, so stages keep reading one
table (``sample_tsv.read_sample_table``) and the user's workbook is never
rewritten. Nextflow reads the table before any process runs and so cannot
take a workbook: run ``EGAP --prepare_only`` on the workbook, then give
Nextflow the derived TSV.

The columns of sheets 1 and 2, their required/format/help text and their
allowed values all come from ``sample_tsv.SAMPLE_COLUMNS``.

Scope of sheet 3, "Sample Settings" (follow-on work, not implemented)
---------------------------------------------------------------------
* One row per SAMPLE_ID of sheet 1, pre-filled; one column per tool setting,
  titled ``{tool}_{setting_name}`` (``filtlong_min_mean_q``,
  ``flye_iterations``), each with a tooltip giving what the setting does, its
  typical range and EGAP's default.
* One registry (planned ``bin/tool_settings.py``) of
  ``Setting(tool, name, default, type or choices, typical range, help, stage)``
  is the single source for three things: what the stage modules run (they
  hard-code their parameters today), these columns and tooltips, and
  ``EGAP.get_pipeline_settings``. It replaces the display-only ``*_SETTINGS``
  dicts in EGAP.py, several of which no longer match the commands in bin/.
* Precedence: built-in default, then a run-level command-line value, then the
  sample's cell. Values are type and range checked at start-up.
* Stages receive their values through the table they already read: this sheet
  is joined into the derived TSV exactly as the metadata sheet is, so no stage
  signature changes.
* The printed settings gain a per-sample level listing only the overrides
  (``EGAP.print_nested_dict`` already prints to any depth).
* Coverage in two tiers: every parameter EGAP sets today first; then each
  tool's remaining options from its ``--help`` at the pinned version, plus a
  ``{tool}_extra_args`` column so any option is reachable before it has a
  column of its own.
* Not editable: paths EGAP computes, thread and memory values, values derived
  from EST_SIZE (expose the coverage multiplier instead), measured insert
  sizes, and flags that downstream parsing depends on.
* Two hazards to design for: read preprocessing writes into a folder shared by
  all samples of a species, so those samples must agree on preprocessing
  settings; and every stage skips when its output exists, so a changed
  setting does nothing on a rerun into the same output folder.
Until then, a "Sample Settings" sheet is only checked against the SAMPLE_ID rule.

Stage:
    Cross-cutting (input handling, before any stage runs)

Author: Ian Bollinger (ian.bollinger@entheome.org)
"""
import os

import pandas as pd

from sample_tsv import METADATA_SHEET, SAMPLE_COLUMNS, SAMPLE_SHEET, SETTINGS_SHEET

WORKBOOK_EXTENSIONS = (".xlsx", ".xlsm")
_EXCEL_MAX_ROW = 1048576


def is_workbook(path):
    """True when *path* names a workbook rather than a TSV."""
    return str(path).lower().endswith(WORKBOOK_EXTENSIONS)


def _clean(sheet):
    """Strip header and SAMPLE_ID whitespace and drop rows that are entirely blank."""
    sheet = sheet.rename(columns=lambda name: str(name).strip())
    sheet = sheet.dropna(how="all")
    if "SAMPLE_ID" in sheet.columns:
        sheet["SAMPLE_ID"] = sheet["SAMPLE_ID"].str.strip()
    return sheet


def _repeated(values):
    """Return the sorted values that occur more than once."""
    return sorted(values[values.duplicated()].unique())


# --------------------------------------------------------------
# Workbook -> the one table the stages read
# --------------------------------------------------------------
def workbook_to_tsv(workbook_path, output_dir):
    """Validate a sample workbook and write the TSV the stages read.

    Parameters
    ----------
    workbook_path : str
        The ``.xlsx`` sample workbook.
    output_dir : str
        Run output folder; created if missing.

    Returns
    -------
    str
        Path of ``<output_dir>/<workbook stem>_samples.tsv``.

    Raises
    ------
    ValueError
        When the first sheet has no usable SAMPLE_ID column, a SAMPLE_ID is
        blank or repeated, a later sheet names a SAMPLE_ID that the first sheet
        does not have, or the metadata sheet repeats a sample-sheet column.
    """
    sheets = pd.read_excel(workbook_path, sheet_name=None, dtype=str)
    names = list(sheets)
    samples = _clean(sheets[names[0]])
    if "SAMPLE_ID" not in samples.columns:
        raise ValueError(f"The first sheet ('{names[0]}') of {workbook_path} has no SAMPLE_ID column.")
    if samples["SAMPLE_ID"].isna().any() or (samples["SAMPLE_ID"] == "").any():
        raise ValueError(f"Sheet '{names[0]}' has a row with a blank SAMPLE_ID.")
    if _repeated(samples["SAMPLE_ID"]):
        raise ValueError(f"Sheet '{names[0]}' repeats SAMPLE_ID: {', '.join(_repeated(samples['SAMPLE_ID']))}")
    known = set(samples["SAMPLE_ID"])

    keyed = {}
    for name in names[1:]:
        sheet = _clean(sheets[name])
        if "SAMPLE_ID" not in sheet.columns:
            continue  # a notes or terms sheet
        sheet = sheet[sheet["SAMPLE_ID"].notna() & (sheet["SAMPLE_ID"] != "")]
        unknown = sorted(set(sheet["SAMPLE_ID"]) - known)
        if unknown:
            raise ValueError(f"Sheet '{name}' has rows for SAMPLE_ID not on sheet '{names[0]}': "
                             f"{', '.join(unknown)}")
        if _repeated(sheet["SAMPLE_ID"]):
            raise ValueError(f"Sheet '{name}' repeats SAMPLE_ID: {', '.join(_repeated(sheet['SAMPLE_ID']))}")
        keyed[name] = sheet

    table = samples
    if METADATA_SHEET in keyed:
        metadata = keyed[METADATA_SHEET]
        clash = sorted((set(metadata.columns) & set(samples.columns)) - {"SAMPLE_ID"})
        if clash:
            raise ValueError(f"Sheet '{METADATA_SHEET}' repeats column(s) of sheet '{names[0]}': "
                             f"{', '.join(clash)}")
        table = samples.merge(metadata, on="SAMPLE_ID", how="left")
    if SETTINGS_SHEET in keyed:
        print(f"NOTE:\tSheet '{SETTINGS_SHEET}' was checked against the sample list but is not "
              f"read yet; per-sample settings are not applied in this version.")

    os.makedirs(output_dir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(str(workbook_path)))[0]
    tsv_path = os.path.join(output_dir, f"{stem}_samples.tsv")
    table.to_csv(tsv_path, sep="\t", index=False)
    return tsv_path


# --------------------------------------------------------------
# Generate a workbook with a tooltip on every column
# --------------------------------------------------------------
def tooltip_text(column):
    """Return the hover text for one column header."""
    lines = [f"Required: {column.required}", f"Format: {column.format}", column.help]
    if column.choices:
        lines.append("Allowed: " + ", ".join(column.choices))
    return "\n".join(lines)


def write_workbook(workbook_path, rows=None):
    """Write a sample workbook: both sheets, tooltips, dropdowns, frozen header.

    Every header cell gets a comment (shown on hover in Excel, LibreOffice and
    as a note in Google Sheets) built from ``sample_tsv.SAMPLE_COLUMNS``.
    Columns with a fixed set of values also get a dropdown.

    Parameters
    ----------
    workbook_path : str
        Where to write the ``.xlsx``.
    rows : pandas.DataFrame, optional
        Existing samples (the raw, collapsed TSV layout) to pre-fill. Columns
        EGAP does not know are kept on the first sheet with a tooltip saying so.

    Returns
    -------
    str
        *workbook_path*.
    """
    from openpyxl import Workbook
    from openpyxl.comments import Comment
    from openpyxl.styles import Font
    from openpyxl.utils import get_column_letter
    from openpyxl.worksheet.datavalidation import DataValidation

    known = {column.name for column in SAMPLE_COLUMNS}
    extra = [name for name in (rows.columns if rows is not None else [])
             if str(name).strip() and name not in known]
    book = Workbook()
    book.remove(book.active)
    for sheet_name in (SAMPLE_SHEET, METADATA_SHEET):
        sheet = book.create_sheet(sheet_name)
        columns = [(c.name, tooltip_text(c), c.choices) for c in SAMPLE_COLUMNS if c.sheet == sheet_name]
        if sheet_name == SAMPLE_SHEET:
            columns += [(name, "Required: no\nFormat: any\nNot an EGAP column. Carried through "
                               "to the stages unchanged.", ()) for name in extra]
        for index, (name, tooltip, choices) in enumerate(columns, start=1):
            header = sheet.cell(row=1, column=index, value=name)
            header.font = Font(bold=True)
            comment = Comment(tooltip, "EGAP")
            comment.width, comment.height = 340, 30 + 18 * (tooltip.count("\n") + len(tooltip) // 55)
            header.comment = comment
            letter = get_column_letter(index)
            sheet.column_dimensions[letter].width = max(16, len(name) + 4)
            if choices:
                allowed = DataValidation(type="list", formula1='"' + ",".join(choices) + '"',
                                         allow_blank=True, showErrorMessage=True,
                                         errorTitle=name[:32],
                                         error=("Use one of: " + ", ".join(choices))[:255])
                sheet.add_data_validation(allowed)
                allowed.add(f"{letter}2:{letter}{_EXCEL_MAX_ROW}")
            if rows is not None and name in rows.columns:
                for row_number, value in enumerate(rows[name], start=2):
                    if pd.notna(value) and str(value).strip() != "":
                        sheet.cell(row=row_number, column=index, value=str(value))
        sheet.freeze_panes = "A2"
    book.save(workbook_path)
    return str(workbook_path)
