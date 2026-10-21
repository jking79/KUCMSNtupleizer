# Style review tool

`style_review.py` audits tracked ntuplizer and skimmer source files without
including backups, generated products, ROOT files, compressed files, or cached
request data. Its maintained-source scope is defined in
`.style-review.json`.

Start with the read-only inventory:

```bash
python3 tools/style_review.py survey
```

Preview the narrow mechanical cleanup as a patch:

```bash
python3 tools/style_review.py diff
```

`diff` exits with status 1 when it finds changes, matching the behavior of the
CI-oriented `check` command.

Check the same policy in automation:

```bash
python3 tools/style_review.py check
```

Apply the byte-safe cleanup only after reviewing the diff:

```bash
python3 tools/style_review.py apply
```

The built-in cleanup normalizes LF line endings, removes ordinary trailing
spaces and tabs, and adds final newlines. It never changes indentation. Review
the patch before applying it because trailing whitespace can be meaningful in
multiline strings or heredocs. Whitespace after a trailing backslash is always
left unchanged and reported for manual review because removing it can change
C/C++, Make, shell, or Python continuation behavior.

## External formatters

Pass `--formatters` to `diff`, `check`, or `apply` to additionally use:

- `clang-format` for C and C++ files;
- Black for Python files;
- `shfmt` for shell scripts.

The formatters operate through standard input and cannot edit repository files
directly. The command fails if a required formatter is missing. Pin formatter
versions in CI before enabling this mode as an enforced repository policy.

Run a narrower review by putting global options before the subcommand:

```bash
python3 tools/style_review.py --path plugins --path KUCMSSkimmer/src survey
```

`apply` refuses to run when tracked files are already modified. The
`--allow-dirty` override exists for deliberate use, but a dedicated formatting
commit from a clean worktree is the recommended baseline workflow.
