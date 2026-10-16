import importlib.util
from pathlib import Path
import tempfile
import unittest


MODULE_PATH = Path(__file__).parents[1] / "style_review.py"
SPEC = importlib.util.spec_from_file_location("style_review", MODULE_PATH)
assert SPEC and SPEC.loader
style_review = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(style_review)


class NormalizeBytesTest(unittest.TestCase):
    def test_normalizes_safe_whitespace_without_changing_indentation(self):
        original = b"\tvalue = 1;  \r\n    next();\r\nlast"
        cleaned, issues = style_review.normalize_bytes(original)

        self.assertEqual(cleaned, b"\tvalue = 1;\n    next();\nlast\n")
        self.assertIn("crlf", issues)
        self.assertIn("trailing-whitespace", issues)
        self.assertIn("missing-final-newline", issues)

    def test_preserves_whitespace_after_continuation_backslash(self):
        original = b"value = one \\\t\n    two\n"
        cleaned, issues = style_review.normalize_bytes(original)

        self.assertEqual(cleaned, original)
        self.assertIn("guarded-backslash-whitespace", issues)
        self.assertNotIn("trailing-whitespace", issues)

    def test_preserves_make_recipe_tab(self):
        original = b"target:\n\tcommand\n"
        cleaned, issues = style_review.normalize_bytes(original)

        self.assertEqual(cleaned, original)
        self.assertEqual(issues, ())

    def test_binary_is_not_modified(self):
        original = b"ROOT\x00data\r\n"
        cleaned, issues = style_review.normalize_bytes(original)

        self.assertEqual(cleaned, original)
        self.assertEqual(issues, ("binary",))


class SelectionTest(unittest.TestCase):
    def test_scope_excludes_backups(self):
        config = {
            "include": ["plugins/**", "KUCMSSkimmer/src/**"],
            "exclude": ["**/backups/**", "**/temp/**", "**/*.bak"],
        }
        paths = [
            "plugins/KUCMSNtupleizer.hh",
            "plugins/backups/KUCMSNtupleizer.hh.bak",
            "KUCMSSkimmer/src/KUCMSAodSkimmer_Photons.cc",
            "KUCMSSkimmer/src/temp/generated.cc",
            "macros/legacy.cpp",
        ]

        selected, excluded = style_review.select_paths(paths, config, [])

        self.assertEqual(
            selected,
            [
                "plugins/KUCMSNtupleizer.hh",
                "KUCMSSkimmer/src/KUCMSAodSkimmer_Photons.cc",
            ],
        )
        self.assertIn("plugins/backups/KUCMSNtupleizer.hh.bak", excluded)
        self.assertIn("KUCMSSkimmer/src/temp/generated.cc", excluded)
        self.assertIn("macros/legacy.cpp", excluded)


class AtomicWriteTest(unittest.TestCase):
    def test_preserves_file_mode(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "script.sh"
            path.write_bytes(b"old\n")
            path.chmod(0o755)

            style_review.write_atomic(path, b"new\n")

            self.assertEqual(path.read_bytes(), b"new\n")
            self.assertEqual(path.stat().st_mode & 0o777, 0o755)


if __name__ == "__main__":
    unittest.main()
