"""Validate a bounded sample-label/file-reference CSV without reading read data."""
from __future__ import annotations

import argparse
import csv
import ctypes
import json
import os
from pathlib import Path, PureWindowsPath
import re
import stat
import sys
import unicodedata


class SampleSheetError(Exception):
    """A declared input failure with a stable category for callers and tests."""

    def __init__(self, category: str, explanation: str):
        super().__init__(explanation)
        self.category = category


def check_local_form(value: str) -> None:
    """Reject network, URI and ambiguous Windows forms before path lookup."""
    # File-reference policy is local and does not expand shell expressions (D-011).
    windows = PureWindowsPath(value)
    if re.match(r"^[A-Za-z][A-Za-z0-9+.-]*://", value):
        raise SampleSheetError("path", "URLs are not supported")
    if value.startswith(("\\\\", "//")) or windows.drive.startswith("\\\\"):
        raise SampleSheetError("path", "UNC/network paths are not supported")
    if windows.drive:
        if not re.fullmatch(r"[A-Za-z]:", windows.drive) or not windows.root:
            raise SampleSheetError("path", "use an ordinary absolute path, not a drive-relative path")
    elif windows.root:
        raise SampleSheetError("path", "rooted paths without a drive are not supported")
    if ":" in value[2 if windows.drive else 0:]:
        raise SampleSheetError("path", "URI or alternate-stream path forms are not supported")


def check_local_drive(path: Path) -> None:
    """Reject mapped network drives as well as explicit UNC path spellings."""
    if os.name == "nt":
        drive_type = ctypes.windll.kernel32.GetDriveTypeW
        drive_type.argtypes = [ctypes.c_wchar_p]
        drive_type.restype = ctypes.c_uint
        if drive_type(str(path.anchor)) == 4:
            raise SampleSheetError("path", "mapped network drives are not supported")


def resolve_local(value: str, base: Path) -> Path:
    check_local_form(value)
    path = Path(value)
    if not path.is_absolute():
        path = base / path
    check_local_drive(path)
    target = path.resolve(strict=True)
    check_local_form(str(target))
    check_local_drive(target)
    return target


def check_value(value: str, record: int) -> None:
    if not value or value != value.strip():
        raise SampleSheetError("row", f"record {record}: values must be nonempty without surrounding whitespace")
    if any(unicodedata.category(character) == "Cc" for character in value):
        raise SampleSheetError("row", f"record {record}: control characters are not allowed")


def count_samples(sheet: Path) -> int:
    """Check the whole table and references; labels are not biological identities."""
    try:
        resolved = resolve_local(str(sheet), Path.cwd())
        handle = resolved.open("r", encoding="utf-8-sig", newline="")
    except SampleSheetError:
        raise
    except (OSError, ValueError) as error:
        raise SampleSheetError("sheet", f"cannot open sample sheet: {error}") from error

    labels: set[str] = set()
    targets: list[Path] = []
    reader = None
    try:
        with handle:
            # Use the pinned standard reader, with explicit text/row grammar (D-010).
            reader = csv.reader(handle, dialect="excel", strict=True)
            if next(reader, None) != ["sample_id", "fastq_path"]:
                raise SampleSheetError("header", "expected exactly sample_id,fastq_path in that order")
            for record, row in enumerate(reader, start=2):
                if len(row) != 2:
                    raise SampleSheetError("row", f"record {record}: expected exactly two values")
                label, value = row
                check_value(label, record)
                check_value(value, record)
                if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", label):
                    raise SampleSheetError("id", f"record {record}: invalid sample ID")
                if label in labels:
                    raise SampleSheetError("id", f"record {record}: duplicate sample ID {label!r}")
                try:
                    target = resolve_local(value, resolved.parent)
                    if not stat.S_ISREG(target.stat().st_mode):
                        raise SampleSheetError("path", f"record {record}: reference is not a regular file")
                    # Metadata acceptance does not inspect FASTQ payloads (D-011).
                    with target.open("rb"):
                        pass
                    if any(target.samefile(previous) for previous in targets):
                        raise SampleSheetError("reuse", f"record {record}: the same underlying file is reused")
                except SampleSheetError as error:
                    if not str(error).startswith("record "):
                        raise SampleSheetError(error.category, f"record {record}: {error}") from error
                    raise
                except FileNotFoundError as error:
                    raise SampleSheetError("path", f"record {record}: referenced file is missing") from error
                except (OSError, ValueError) as error:
                    raise SampleSheetError("access", f"record {record}: cannot access or identify reference: {error}") from error
                labels.add(label)
                targets.append(target)
    except UnicodeError as error:
        raise SampleSheetError("encoding", "sample sheet is not valid UTF-8") from error
    except csv.Error as error:
        line = reader.line_num if reader is not None else 0
        raise SampleSheetError("csv", f"CSV grammar error near physical line {line}: {error}") from error
    except OSError as error:
        raise SampleSheetError("sheet", f"cannot read sample sheet: {error}") from error
    if not targets:
        raise SampleSheetError("row", "sample sheet must contain at least one data row")
    return len(targets)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, color=False)
    parser.add_argument("sheet", type=Path, help="UTF-8 CSV with sample_id,fastq_path columns")
    args = parser.parse_args(argv)
    try:
        count = count_samples(args.sheet)
    except SampleSheetError as error:
        print(f"{error.category}: {error}", file=sys.stderr)
        return 2
    # Emit only after complete acceptance; never a partial biological count (D-012).
    print(json.dumps({"samples": count}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
