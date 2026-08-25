# tests/gates/baseline — G0.1 identity instrument

The off-path promise: `--smooth` absent is **byte-identical to shipped
1.0.0**. This directory is the machine that proves it. It is not a drawer
of golden `.step` files.

**Why `187ead0`.** That commit is released stl2step 1.0.0 (`stl2step
--version` → `stl2step 1.0.0`), the last engine before any smooth work.
G0.1 compares live output of a P0-built binary of that commit against the
current CLI. Never commit the binary or the STEP bytes; rebuild it.

**Rebuild.** From the repo root (or here):

```
tests/gates/baseline/build_baseline.sh
```

Idempotent. `--force` wipes and starts over. Prints `BASELINE_BIN=...`
and `--version`. Configure/build type follow `<repo>/build` when that
cache exists, otherwise Unix Makefiles + Release.

**What is compared (and what is stripped).**

| Surface | Canonical form | Stripped, and why |
|---|---|---|
| STEP | `DATA; … ENDSEC;` only | Whole HEADER, because `FILE_NAME` field 2 is a wall-clock timestamp. Other HEADER fields (description, name, author, organisation, preprocessor version, originating system, authorisation, `FILE_SCHEMA`) are *not* rewritten individually — they live in HEADER, so DATA-only compare drops them. Do **not** strip anything in DATA; `#N` entity ids are compared verbatim. |
| RESULT | ordered key list + values | `seconds` (wall clock) and `input`/`output` values (resolved paths). Keys themselves, including those two, stay in the ordered list. **No `smooth*` / `facesAfterSmooth` key may appear.** |

```
python3 tests/gates/baseline/canonicalize.py step  a.step b.step
python3 tests/gates/baseline/canonicalize.py result a.txt b.txt
# exit 0 identical, 1 different (unified diff of canonical forms), 2 IO/usage
```

p0-gates wires this with `--baseline tests/gates/baseline`: build the
187ead0 CLI via `build_baseline.sh`, twin-run it against the current
binary on the same STL, then call `canonicalize.py`. There are no
checked-in `step/` or `result/` goldens.

**Warning.** Loosening a strip rule to make a gate green is a defect, not
a fix. If current-vs-baseline fails, report the field that moved; do not
teach the canonicalizer to ignore it.
