"""P2 tests against manual quality profiles, not runner-generated expectations."""
from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
import zipfile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import report_read_quality as reporter

FIXTURES = ROOT / "tests/fixtures/read_quality"
CLI = ROOT / "scripts/report_read_quality.py"


def parse_report(archive_path: Path) -> dict:
    """Read the existing FastQC data format independently of the runner."""
    with zipfile.ZipFile(archive_path) as archive:
        name = next(n for n in archive.namelist() if n.endswith("/fastqc_data.txt"))
        lines = archive.read(name).decode("utf-8").splitlines()
    result = {"version": lines[0].split("\t")[1], "modules": {}, "basic": {}, "quality": []}
    module = None
    for line in lines[1:]:
        if line == ">>END_MODULE":
            module = None
        elif line.startswith(">>"):
            module, flag = line[2:].split("\t")
            result["modules"][module] = flag
        elif line.startswith("#") or not line:
            continue
        elif module == "Basic Statistics":
            key, value = line.split("\t", 1)
            if key != "Filename":
                result["basic"][key] = value
        elif module == "Per base sequence quality":
            values = line.split("\t")
            result["quality"].append([int(values[0]), *map(float, values[1:])])
    return result


class ReadQualityTests(unittest.TestCase):
    def test_fixture_bytes_and_manual_profiles(self):
        # Exact fixture contents/oracles precede actual evaluation (D-008).
        hashes = json.loads((FIXTURES / "sha256.json").read_text())
        for name, size, quality in (("control", 9797, b"?" * 40), ("declining", 9999, b"?" * 20 + b"+" * 20)):
            with self.subTest(name=name):
                path = FIXTURES / (name + ".fastq")
                data = path.read_bytes()
                self.assertEqual(len(data), size)
                self.assertEqual(hashlib.sha256(data).hexdigest(), hashes[path.name])
                self.assertNotIn(b"\r", data)
                self.assertTrue(data.endswith(b"\n"))
                lines = data.splitlines()
                self.assertEqual(len(lines), 404)
                for index in range(101):
                    self.assertEqual(lines[4 * index], f"@{name}_{index + 1:03d}".encode("ascii"))
                    self.assertEqual(lines[4 * index + 1], b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT")
                    self.assertEqual(lines[4 * index + 2], b"+")
                    self.assertEqual(lines[4 * index + 3], quality)
                self.assertEqual(reporter.count_fastq(path), 101)

    def test_real_input_setup_and_output_refusals(self):
        # Input/setup rejection must precede analysis and preserve old outputs (D-006).
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            malformed = base / "bad.fastq"
            malformed.write_bytes(b"@bad\nACGT\n+\nIII\n")
            existing = base / "existing"
            existing.mkdir()
            marker = existing / "keep.txt"
            marker.write_text("unchanged")
            cases = [
                (base / "missing.fastq", base / "missing-out", CLI, Path(sys.executable)),
                (malformed, base / "bad-out", CLI, Path(sys.executable)),
                (FIXTURES / "control.fastq", base / "no-java", CLI, base / "missing-java"),
                (FIXTURES / "control.fastq", base / "no-fastqc", base / "missing-fastqc", Path(sys.executable)),
                (FIXTURES / "control.fastq", existing, CLI, Path(sys.executable)),
            ]
            for source, output, fastqc, java in cases:
                with self.subTest(output=output.name):
                    result = subprocess.run([sys.executable, str(CLI), str(source), str(output),
                                             "--fastqc", str(fastqc), "--java", str(java)], capture_output=True)
                    self.assertEqual(result.returncode, 2, result.stderr)
                    self.assertEqual(result.stdout, b"")
                    self.assertTrue(result.stderr)
                    if output != existing:
                        self.assertFalse(output.exists())
            self.assertEqual(marker.read_text(), "unchanged")

    def test_injected_tool_and_incomplete_output_failures(self):
        # These explicitly injected failures test handling, not real FastQC errors.
        for tool_code in (7, 0):
            with self.subTest(tool_code=tool_code), tempfile.TemporaryDirectory() as temporary:
                base = Path(temporary)
                launcher = base / "fastqc"
                launcher.write_text("fixture launcher placeholder")
                (base / "Configuration").mkdir()
                (base / "Configuration/limits.txt").write_text("fixture config")
                output = base / "output"

                def fake_command(command, directory, label):
                    (directory / f"{label}.stdout.txt").write_text("injected")
                    (directory / f"{label}.stderr.txt").write_text("injected")
                    text = {"java-version": "Temurin 21.0.13.1+2", "fastqc-version": "FastQC v0.13.0"}.get(label, "")
                    if label == "fastqc":
                        (directory / "partial.txt").write_text("retained partial evidence")
                    return {"command": command, "exit_code": tool_code if label == "fastqc" else 0,
                            "seconds": 0, "stdout": text, "stderr": ""}

                with patch.object(reporter, "run_command", side_effect=fake_command), redirect_stderr(io.StringIO()), redirect_stdout(io.StringIO()):
                    code = reporter.main([str(FIXTURES / "control.fastq"), str(output),
                                          "--fastqc", str(launcher), "--java", sys.executable])
                self.assertEqual(code, 1)
                metadata = json.loads((output / "run.json").read_text())
                self.assertEqual(metadata["status"], "failed")
                self.assertEqual(metadata["failure_phase"], "analysis")
                self.assertTrue((output / "fastqc.stderr.txt").exists())
                self.assertTrue((output / "partial.txt").exists())

    def test_cli_help_and_missing_arguments(self):
        help_result = subprocess.run([sys.executable, str(CLI), "--help"], capture_output=True)
        missing = subprocess.run([sys.executable, str(CLI)], capture_output=True)
        self.assertEqual(help_result.returncode, 0)
        self.assertIn(b"FastQC", help_result.stdout)
        self.assertEqual(missing.returncode, 2)
        self.assertEqual(missing.stdout, b"")

    @unittest.skipUnless(os.environ.get("P2_FASTQC") and os.environ.get("P2_JAVA"), "set explicit pinned FastQC/Java paths for real integration")
    def test_real_fastqc_profiles_and_execution_status(self):
        retained = os.environ.get("P2_REPORT_EVIDENCE_DIR")
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(retained) if retained else Path(temporary)
            base.mkdir(parents=True, exist_ok=True)
            for name in ("control", "declining"):
                with self.subTest(name=name):
                    source = FIXTURES / (name + ".fastq")
                    before = reporter.sha256(source)
                    output = base / name
                    result = subprocess.run([sys.executable, str(CLI), str(source), str(output),
                                             "--fastqc", os.environ["P2_FASTQC"], "--java", os.environ["P2_JAVA"]], capture_output=True)
                    self.assertEqual(result.returncode, 0, result.stderr.decode(errors="replace"))
                    self.assertEqual(reporter.sha256(source), before)
                    metadata = json.loads((output / "run.json").read_text())
                    self.assertEqual(metadata["status"], "complete")
                    self.assertTrue(list(output.glob("*_fastqc.html")))
                    report = parse_report(next(output.glob("*_fastqc.zip")))
                    self.assertEqual(report["version"], "0.13.0")
                    self.assertEqual(report["basic"]["Total Sequences"], "101")
                    self.assertEqual(report["basic"]["Sequence length"], "40")
                    self.assertEqual(report["basic"]["Encoding"], "Sanger / Illumina 1.9")
                    self.assertEqual(len(report["quality"]), 40)
                    for position, row in enumerate(report["quality"], 1):
                        expected = 10 if name == "declining" and position > 20 else 30
                        self.assertEqual(row, [position] + [expected] * 6)
                    self.assertEqual(report["modules"]["Per base sequence quality"], "pass" if name == "control" else "fail")
                    if retained:
                        (base / (name + "-parsed.json")).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    unittest.main()
