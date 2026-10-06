"""Regression checks for certificate verification and its published transcript."""

import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import verify


ROOT = Path(__file__).resolve().parents[1]


class VerifierTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.parameters = json.loads((ROOT / "data/182-vertex-code.json").read_text())

    def test_rebuilds_published_summary_and_every_transcript_entry(self):
        with tempfile.TemporaryDirectory() as directory:
            trace = Path(directory) / "transcript.json"
            result = verify.verify(self.parameters, trace)
            # JSON turns histogram keys into strings; compare in that format.
            self.assertEqual(
                json.loads(json.dumps(result)),
                json.loads((ROOT / "data/verification.json").read_text()),
            )
            self.assertEqual(
                json.loads(trace.read_text()),
                json.loads((ROOT / "data/reference-transcript.json").read_text()),
            )

    def test_rejects_malformed_finite_parameters(self):
        mutants = []
        for name in ("order", "backward", "C", "J"):
            duplicate = copy.deepcopy(self.parameters)
            duplicate[name][0] = duplicate[name][1]
            mutants.append((name + " duplicate", duplicate))
            for value in (True, 1.0, "1", -1, 910):
                invalid = copy.deepcopy(self.parameters)
                invalid[name][0] = value
                mutants.append((name + " invalid " + repr(value), invalid))
        extra = copy.deepcopy(self.parameters)
        extra["ignored"] = []
        mutants.append(("extra key", extra))
        missing = copy.deepcopy(self.parameters)
        del missing["J"]
        mutants.append(("missing key", missing))
        overlap = copy.deepcopy(self.parameters)
        overlap["J"][0] = overlap["C"][0]
        mutants.append(("C/J overlap", overlap))
        for description, data in mutants:
            with self.subTest(description=description):
                with self.assertRaises(verify.VerificationError):
                    verify.verify(data)

    def test_rejects_duplicate_json_keys(self):
        with self.assertRaisesRegex(verify.VerificationError, "Duplicate JSON key"):
            json.loads('{"order": [], "order": []}',
                       object_pairs_hook=verify.unique_object)

    def test_checks_survive_optimization_and_unrelated_working_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            success = subprocess.run(
                [sys.executable, "-O", str(ROOT / "verify.py")],
                cwd=directory, capture_output=True, text=True, check=False,
            )
            self.assertEqual(success.returncode, 0, success.stderr)
            self.assertTrue(json.loads(success.stdout)["all_decoders_verified"])
            invalid = copy.deepcopy(self.parameters)
            invalid["order"][0] = invalid["order"][1]
            path = Path(directory) / "invalid.json"
            path.write_text(json.dumps(invalid))
            failure = subprocess.run(
                [sys.executable, "-O", str(ROOT / "verify.py"), str(path)],
                cwd=directory, capture_output=True, text=True, check=False,
            )
            self.assertNotEqual(failure.returncode, 0)
            self.assertIn("Verification failed", failure.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
