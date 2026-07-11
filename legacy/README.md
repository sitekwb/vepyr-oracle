# legacy/

READ-ONLY vendored copies of the original ad-hoc validation harness that
lived untracked at `~/vepyr/work/` on the HPC cluster. They are kept here
purely for reference/diffing while the pipeline is decomposed into the
`oracle/` package proper (see later tasks). Do not import from or extend
these files — port their logic into `oracle/*.py` instead.

Notes:
- `validate.py` does `import vepyr` at module top, so it cannot be imported
  or unit-tested outside the cluster environment as-is.
- `diff_test.py` was dropped as stale: it does `from validate import diff`,
  but `validate.py` only defines `diff_stream` (no `diff` function) — the
  import would fail.
