"""P1 contract tests with predeclared manual outcomes, not parser-derived answers."""
from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
CLI = ROOT / "scripts/validate_fastq.py"
FIXTURES = ROOT / "tests/fixtures/fastq"
CASES = json.loads((FIXTURES / "expectations.json").read_text(encoding="utf-8"))
spec = importlib.util.spec_from_file_location("validate_fastq", CLI)
validator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validator)


def retain_result(name: str, code: int, stdout: bytes, stderr: bytes) -> None:
    """Optionally retain observed bytes for the two-environment comparison."""
    # Keep actual per-case results separate from predeclared expectations (D-005).
    target = os.environ.get("P1_CASE_EVIDENCE_DIR")
    if not target:
        return
    directory = Path(target)
    directory.mkdir(parents=True, exist_ok=True)
    if (directory / f"{name}.json").exists():
        raise AssertionError(f"refusing to overwrite prior case evidence: {name}")
    (directory / f"{name}.stdout.bin").write_bytes(stdout)
    (directory / f"{name}.stderr.bin").write_bytes(stderr)
    category = stderr.decode("utf-8", errors="replace").split(":", 1)[0] if stderr else None
    (directory / f"{name}.json").write_text(
        json.dumps({"exit": code, "category": category}, indent=2) + "\n", encoding="utf-8"
    )


class ValidateFASTQTests(unittest.TestCase):
    def run_cli(self, path: Path):
        return subprocess.run([sys.executable, str(CLI), str(path)], capture_output=True, timeout=15)

    def test_predeclared_cli_cases_and_immutable_inputs(self):
        # The outcomes instantiate the previously written fixture contract (D-005).
        for case in CASES:
            with self.subTest(case=case["id"]):
                path = FIXTURES / case["file"]
                before = hashlib.sha256(path.read_bytes()).hexdigest()
                result = self.run_cli(path)
                retain_result(case["id"], result.returncode, result.stdout, result.stderr)
                self.assertEqual(result.returncode, case["exit"], result.stderr)
                self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), before)
                if case["exit"] == 0:
                    self.assertEqual(json.loads(result.stdout), {"reads": case["reads"]})
                    self.assertTrue(result.stdout.endswith(b"\n"))
                    self.assertEqual(result.stderr, b"")
                else:
                    # A valid prefix cannot leak a partial success count (D-002).
                    self.assertEqual(result.stdout, b"")
                    self.assertTrue(result.stderr.startswith(case["category"].encode("ascii") + b":"))
                    self.assertNotIn(b"Traceback", result.stderr)

    def test_missing_path_is_a_controlled_input_access_error(self):
        with tempfile.TemporaryDirectory(prefix="p1 missing ") as temporary:
            path = Path(temporary) / "absent.fastq"
            self.assertFalse(path.exists())
            result = self.run_cli(path)
            retain_result("N01", result.returncode, result.stdout, result.stderr)
        self.assertEqual(result.returncode, 2)
        self.assertEqual(result.stdout, b"")
        self.assertTrue(result.stderr.startswith(b"input access:"))

    def test_unreadable_input_simulation_is_reported(self):
        # Simulated opening failure is approved; this is not an actual ACL test.
        stdout, stderr = io.StringIO(), io.StringIO()
        with patch.object(Path, "open", side_effect=PermissionError(13, "permission denied")):
            with redirect_stdout(stdout), redirect_stderr(stderr):
                code = validator.main([str(FIXTURES / "P01.fastq")])
        retain_result("N02-focused-simulation", code, stdout.getvalue().encode(), stderr.getvalue().encode())
        self.assertEqual(code, 2)
        self.assertEqual(stdout.getvalue(), "")
        self.assertTrue(stderr.getvalue().startswith("input access:"))

    def test_help_is_successful_usage_without_a_count(self):
        result = subprocess.run([sys.executable, str(CLI), "--help"], capture_output=True, timeout=15)
        retain_result("help", result.returncode, result.stdout, result.stderr)
        self.assertEqual(result.returncode, 0)
        self.assertIn(b"usage:", result.stdout)
        self.assertNotIn(b'"reads"', result.stdout)
        self.assertEqual(result.stderr, b"")

    def test_missing_argument_is_a_usage_error(self):
        result = subprocess.run([sys.executable, str(CLI)], capture_output=True, timeout=15)
        retain_result("missing-argument", result.returncode, result.stdout, result.stderr)
        self.assertEqual(result.returncode, 2)
        self.assertEqual(result.stdout, b"")
        self.assertIn(b"usage:", result.stderr)


if __name__ == "__main__":
    unittest.main()
