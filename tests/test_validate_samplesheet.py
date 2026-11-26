"""P3 manual contract oracles; retained observations are not expected values."""
from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "tests/fixtures/samplesheet"
CLI = ROOT / "scripts/validate_samplesheet.py"
sys.path.insert(0, str(ROOT / "scripts"))
import validate_samplesheet as validator


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SampleSheetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Manual inputs/oracles were frozen before implementation evaluation (D-013).
        cls.temporary = None
        if os.environ.get("P3_FIXTURE_ROOT"):
            cls.base = Path(os.environ["P3_FIXTURE_ROOT"]).resolve()
        else:
            parent = ROOT / ".cache/p3-tests"
            parent.mkdir(parents=True, exist_ok=True)
            cls.temporary = tempfile.TemporaryDirectory(dir=parent)
            cls.base = Path(cls.temporary.name)
            shutil.copytree(SOURCE, cls.base, dirs_exist_ok=True)
            absolute = str((cls.base / "files/a.fastq").resolve()).encode("utf-8")
            case = cls.base / "P3-A07.csv"
            case.write_bytes(case.read_bytes().replace(b"{ABS_A}", absolute))
            os.link(cls.base / "files/a.fastq", cls.base / "files/hardlink.fastq")
            (cls.base / "elsewhere").mkdir()
            (cls.base / "shadow/files").mkdir(parents=True)
            (cls.base / "shadow/files/a.fastq").write_bytes(b"shadow\n")
        cls.specification = json.loads((SOURCE / "expectations.json").read_text(encoding="utf-8"))
        cls.before = {p.relative_to(cls.base).as_posix(): digest(p) for p in cls.base.rglob("*") if p.is_file()}
        cls.observations = {}
        cls.evidence = Path(os.environ["P3_EVIDENCE_DIR"]) if os.environ.get("P3_EVIDENCE_DIR") else None
        if cls.evidence:
            cls.evidence.mkdir(parents=True, exist_ok=True)
            if (cls.evidence / "observations.json").exists():
                raise RuntimeError("refusing to replace retained P3 observations")

    @classmethod
    def tearDownClass(cls):
        after = {p.relative_to(cls.base).as_posix(): digest(p) for p in cls.base.rglob("*") if p.is_file()}
        if after != cls.before:
            raise AssertionError("sample-sheet/reference fixture bytes changed during tests")
        if cls.evidence:
            (cls.evidence / "observations.json").write_text(json.dumps(cls.observations, indent=2) + "\n", encoding="utf-8")
            (cls.evidence / "input-sha256.json").write_text(json.dumps(after, indent=2) + "\n", encoding="utf-8")
        if cls.temporary:
            cls.temporary.cleanup()

    def retain(self, identifier, result, *, simulated=False):
        category = result.stderr.decode("utf-8", errors="replace").split(":", 1)[0] if result.returncode else None
        observation = {"exit_code": result.returncode, "stdout_hex": result.stdout.hex(),
                       "category": category, "simulated": simulated}
        self.observations[identifier] = observation
        if self.evidence:
            target = self.evidence / identifier
            target.mkdir(exist_ok=False)
            (target / "stdout.bin").write_bytes(result.stdout)
            (target / "stderr.bin").write_bytes(result.stderr)
            (target / "result.json").write_text(json.dumps(observation, indent=2) + "\n", encoding="utf-8")

    def invoke(self, sheet, *, cwd=None):
        return subprocess.run([sys.executable, str(CLI), str(sheet)], cwd=cwd, capture_output=True)

    def assert_result(self, result, count=None, category=None):
        if count is not None:
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stderr, b"")
            self.assertEqual(json.loads(result.stdout), {"samples": count})
            self.assertIs(type(json.loads(result.stdout)["samples"]), int)
            self.assertTrue(result.stdout.endswith(b"\n"))
        else:
            self.assertEqual(result.returncode, 2, result.stderr)
            self.assertEqual(result.stdout, b"")
            self.assertTrue(result.stderr.startswith((category + ":").encode()), result.stderr)

    def test_materialized_fixture_bytes_and_file_identity(self):
        hashes = json.loads((SOURCE / "sha256.json").read_text(encoding="utf-8"))
        for name, expected in hashes.items():
            self.assertEqual(digest(SOURCE / name), expected)
            expected_bytes = (SOURCE / name).read_bytes()
            if name == "P3-A07.csv":
                expected_bytes = expected_bytes.replace(b"{ABS_A}", str((self.base / "files/a.fastq").resolve()).encode("utf-8"))
            self.assertEqual((self.base / name).read_bytes(), expected_bytes)
        self.assertTrue((self.base / "files/a.fastq").samefile(self.base / "files/hardlink.fastq"))
        self.assertFalse((self.base / "files/a.fastq").samefile(self.base / "files/copy.fastq"))
        self.assertEqual((self.base / "files/a.fastq").read_bytes(), b"A\n")
        self.assertEqual((self.base / "files/copy.fastq").read_bytes(), b"A\n")

    def test_real_cli_predeclared_table_and_reference_cases(self):
        for case in self.specification["cases"]:
            with self.subTest(case=case["id"]):
                result = self.invoke(self.base / case["file"])
                self.retain(case["id"], result)
                self.assert_result(result, case["count"], case["category"])

    def test_real_working_directory_independence(self):
        for name in ("elsewhere", "shadow"):
            with self.subTest(cwd=name):
                result = self.invoke(self.base / "P3-A02.csv", cwd=self.base / name)
                self.retain("P3-CWD01-" + name, result)
                self.assert_result(result, count=2)

    def test_missing_sheet(self):
        result = self.invoke(self.base / "absent.csv")
        self.retain("P3-IO01", result)
        self.assert_result(result, category="sheet")

    def test_explicitly_injected_access_failures(self):
        original_open, original_stat = Path.open, Path.stat
        for identifier in ("P3-IO02", "P3-IO03", "P3-IO04", "P3-IO05"):
            with self.subTest(case=identifier):
                def open_failure(path, *args, **kwargs):
                    if identifier == "P3-IO02" and path.name == "P3-A02.csv":
                        raise PermissionError("simulated sheet access denial")
                    if identifier == "P3-IO03" and path.name == "a.fastq":
                        raise PermissionError("simulated target access denial")
                    return original_open(path, *args, **kwargs)

                def stat_failure(path, *args, **kwargs):
                    if identifier == "P3-IO04" and path.name == "a.fastq":
                        raise PermissionError("simulated target stat denial")
                    return original_stat(path, *args, **kwargs)

                def identity_failure(path, other):
                    raise OSError("simulated identity lookup failure")

                stdout, stderr = io.StringIO(), io.StringIO()
                identity_patch = patch.object(Path, "samefile", identity_failure) if identifier == "P3-IO05" else patch.object(Path, "stat", stat_failure)
                with patch.object(Path, "open", open_failure), identity_patch, redirect_stdout(stdout), redirect_stderr(stderr):
                    code = validator.main([str(self.base / "P3-A02.csv")])
                result = subprocess.CompletedProcess([], code, stdout.getvalue().encode(), stderr.getvalue().encode())
                self.retain(identifier, result, simulated=True)
                self.assert_result(result, category="sheet" if identifier == "P3-IO02" else "access")

    def test_target_payloads_are_never_read(self):
        original_open = Path.open

        class NoPayloadReads:
            def __init__(self, handle):
                self.handle = handle

            def __enter__(self):
                return self

            def __exit__(self, *args):
                self.handle.close()

            def read(self, *args):
                raise AssertionError("validator attempted to read referenced payload")

            def __iter__(self):
                raise AssertionError("validator attempted to iterate referenced payload")

        def guarded_open(path, mode="r", *args, **kwargs):
            handle = original_open(path, mode, *args, **kwargs)
            return NoPayloadReads(handle) if mode == "rb" else handle

        stdout, stderr = io.StringIO(), io.StringIO()
        with patch.object(Path, "open", guarded_open), redirect_stdout(stdout), redirect_stderr(stderr):
            code = validator.main([str(self.base / "P3-A02.csv")])
        result = subprocess.CompletedProcess([], code, stdout.getvalue().encode(), stderr.getvalue().encode())
        self.retain("P3-CONTENT01", result)
        self.assert_result(result, count=2)

    def test_cli_help_and_missing_argument(self):
        for identifier, arguments, code in [("P3-CLI01", ["--help"], 0), ("P3-CLI02", [], 2)]:
            result = subprocess.run([sys.executable, str(CLI), *arguments], capture_output=True)
            self.retain(identifier, result)
            self.assertEqual(result.returncode, code)
            if code:
                self.assertEqual(result.stdout, b"")
            else:
                self.assertIn(b"sample_id", result.stdout)


if __name__ == "__main__":
    unittest.main()
