"""Validate P1's declared FASTQ contract and report a complete-file read count."""
from __future__ import annotations

import argparse
import io
import json
from pathlib import Path
import sys


class FASTQValidationError(Exception):
    """An input problem with a stable category for callers and diagnostics."""

    def __init__(self, category: str, explanation: str):
        super().__init__(explanation)
        self.category = category


def count_fastq(path: Path) -> int:
    """Check every record without modifying input or retaining previous records."""
    # Use the approved installed reader; no home-grown FASTQ parsing (D-003).
    from Bio import SeqIO

    count = 0
    try:
        with path.open("rb") as raw:
            # Recognize common compressed signatures rather than file suffixes.
            # Compressed input is outside the accepted contract (D-001).
            prefix = raw.read(6)
            signatures = (b"\x1f\x8b", b"BZh", b"\xfd7zXZ\x00", b"PK\x03\x04", b"\x28\xb5\x2f\xfd")
            if any(prefix.startswith(signature) for signature in signatures):
                raise FASTQValidationError("unsupported compression", "provide an uncompressed FASTQ file")
            raw.seek(0)
            with io.TextIOWrapper(raw, encoding="ascii", newline=None) as handle:
                # Reader grammar handles wrapping, captions and qualities (D-001).
                records = iter(SeqIO.parse(handle, "fastq"))
                while True:
                    try:
                        record = next(records)
                    except StopIteration:
                        break
                    except IndexError as error:
                        # Biopython 1.88 indexes the first token of the title;
                        # an empty title must be a controlled input error.
                        raise FASTQValidationError("invalid FASTQ", f"record {count + 1} has an empty identifier") from error
                    if not record.id.strip():
                        raise FASTQValidationError("invalid FASTQ", f"record {count + 1} has an empty identifier")
                    sequence = str(record.seq)
                    if not sequence:
                        raise FASTQValidationError("invalid FASTQ", f"record {count + 1} has an empty sequence")
                    if any(base not in "ACGTNacgtn" for base in sequence):
                        raise FASTQValidationError("invalid FASTQ", f"record {count + 1} has a sequence symbol outside A/C/G/T/N")
                    # Repeated identifiers count separately; no quality cutoff.
                    count += 1
    except UnicodeDecodeError as error:
        raise FASTQValidationError("non-ASCII input", "FASTQ input must contain only ASCII bytes") from error
    except OSError as error:
        raise FASTQValidationError("input access", f"cannot read the input file: {error.strerror or error}") from error
    except ValueError as error:
        raise FASTQValidationError("invalid FASTQ", str(error)) from error
    if count == 0:
        raise FASTQValidationError("empty input", "at least one FASTQ record is required")
    return count


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Check one uncompressed ASCII Phred+33 FASTQ file and count its records.",
        color=False,
    )
    parser.add_argument("fastq", type=Path, help="path to one uncompressed FASTQ file")
    args = parser.parse_args(argv)
    # Emit no count until the complete input has passed; failures stay visible (D-002).
    try:
        count = count_fastq(args.fastq)
    except FASTQValidationError as error:
        print(f"{error.category}: {error}", file=sys.stderr)
        return 2
    print(json.dumps({"reads": count}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
