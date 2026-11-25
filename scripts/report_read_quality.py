"""Validate one FASTQ and retain a pinned FastQC report without filtering reads."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time
import zipfile

from validate_fastq import FASTQValidationError, count_fastq


def sha256(path: Path) -> str:
    """Hash the exact input, configuration or artifact bytes for provenance."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_command(command: list[str], output: Path, label: str) -> dict:
    """Run without a shell and retain both streams, including failed attempts."""
    start = time.perf_counter()
    result = subprocess.run(command, capture_output=True, check=False)
    (output / f"{label}.stdout.txt").write_bytes(result.stdout)
    (output / f"{label}.stderr.txt").write_bytes(result.stderr)
    return {
        "command": command, "exit_code": result.returncode,
        "seconds": time.perf_counter() - start,
        "stdout": result.stdout.decode("utf-8", errors="replace"),
        "stderr": result.stderr.decode("utf-8", errors="replace"),
    }


def report_artifacts(output: Path) -> list[Path]:
    """Require the HTML and intact ZIP with nonempty numeric data and summary."""
    htmls = list(output.glob("*_fastqc.html"))
    archives = list(output.glob("*_fastqc.zip"))
    if len(htmls) != 1 or len(archives) != 1 or not htmls[0].stat().st_size:
        raise ValueError("required HTML/ZIP report artifacts are missing or ambiguous")
    if htmls[0].stem != archives[0].stem:
        raise ValueError("HTML/ZIP report names do not match")
    with zipfile.ZipFile(archives[0]) as archive:
        if archive.testzip() is not None:
            raise ValueError("report ZIP failed its integrity check")
        for name in ("fastqc_data.txt", "summary.txt"):
            matches = [p for p in archive.namelist() if p.endswith("/" + name)]
            if len(matches) != 1 or not archive.read(matches[0]).strip():
                raise ValueError(f"required ZIP content is missing or empty: {name}")
    return [htmls[0], archives[0]]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Check one Phred+33 FASTQ and generate a pinned FastQC report.",
        color=False,
    )
    parser.add_argument("fastq", type=Path)
    parser.add_argument("output", type=Path, help="new output directory; existing paths are refused")
    parser.add_argument("--fastqc", type=Path, required=True, help="FastQC 0.13.0 Python launcher file")
    parser.add_argument("--java", type=Path, required=True, help="portable Temurin java.exe path")
    args = parser.parse_args(argv)
    source, output, launcher, java = [p.resolve() for p in (args.fastq, args.output, args.fastqc, args.java)]
    # Reuse the complete P1 check before downstream work; never repair reads (D-006).
    try:
        if args.output.exists() or args.output.is_symlink():
            raise ValueError("output path already exists; choose a new directory")
        count_fastq(source)
        if not launcher.is_file() or not java.is_file():
            raise ValueError("FastQC launcher and Java executable must be existing files")
        original_hash = sha256(source)
        output.mkdir(parents=True, exist_ok=False)
    except (FASTQValidationError, OSError, ValueError) as error:
        print(f"input/setup: {error}", file=sys.stderr)
        return 2

    metadata = {
        "status": "running", "input": str(source), "input_sha256": original_hash,
        "python": sys.version, "fastqc_launcher": str(launcher), "java": str(java),
        "commands": [],
    }
    code = 1
    phase = "setup"
    try:
        metadata["launcher_sha256"] = sha256(launcher)
        # Fixed artifacts/options and explicit runtime avoid global installation (D-007).
        for label, command in (
            ("java-version", [str(java), "-version"]),
            ("fastqc-version", [sys.executable, str(launcher), "--java", str(java), "--version"]),
        ):
            observation = run_command(command, output, label)
            metadata["commands"].append(observation)
            if observation["exit_code"]:
                raise ValueError(f"{label} failed; see retained logs")
        java_version = metadata["commands"][0]["stdout"] + metadata["commands"][0]["stderr"]
        fastqc_version = metadata["commands"][1]["stdout"].strip()
        if "21.0.13.1+2" not in java_version or "Temurin" not in java_version:
            raise ValueError("Java version does not match approved Temurin 21.0.13.1+2")
        if fastqc_version != "FastQC v0.13.0":
            raise ValueError("FastQC version does not match approved 0.13.0")
        configuration = launcher.parent / "Configuration"
        metadata["configuration_sha256"] = {
            p.name: sha256(p) for p in sorted(configuration.iterdir()) if p.is_file()
        }
        phase = "analysis"
        command = [sys.executable, str(launcher), "--java", str(java), "--format", "fastq",
                   "--nogroup", "--threads", "1", "--memory", "256", "--outdir", str(output), str(source)]
        observation = run_command(command, output, "fastqc")
        metadata["commands"].append(observation)
        if observation["exit_code"]:
            raise ValueError("FastQC execution failed; see retained logs")
        artifacts = report_artifacts(output)
        if sha256(source) != original_hash:
            raise ValueError("input bytes changed during reporting")
        metadata["artifacts"] = {p.name: sha256(p) for p in artifacts}
        # Module WARN/FAIL is report content, not execution failure (D-006).
        metadata["status"] = "complete"
        code = 0
    except (OSError, ValueError, zipfile.BadZipFile) as error:
        code = 2 if phase == "setup" else 1
        metadata.update(status="failed", failure_phase=phase, error=str(error))
        print(f"{phase}: {error}", file=sys.stderr)
    finally:
        (output / "run.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    if code == 0:
        print(f"Report complete: {output}")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
