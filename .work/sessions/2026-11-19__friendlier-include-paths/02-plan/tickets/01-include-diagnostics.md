# Ticket 01: include-diagnostics

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] New `radiopyo/parser/exceptions.py` (banner layout) with `ConfigError(RadioPyoError, ValueError)` carrying an `origin` (file path or `<dict>`, plus include chain), `IncludeNotFoundError(ConfigError, FileNotFoundError)`, `ConfigSyntaxError(ConfigError)`; exported from `radiopyo.parser` and `radiopyo`.
- [ ] `includes = ["a.toml", ...]` (top-level list, or a single string) is accepted; the parser is picked from each entry's suffix.
- [ ] An unknown extension key in the table form, or an entry with an unknown suffix in the list form, raises `ConfigError` naming the file and the accepted forms (no more skipped include).
- [ ] A missing include raises `IncludeNotFoundError` whose message names the including file, the include chain (`main.toml -> a.toml`), the entry as written, every location tried, and close matches (difflib) among `bundled_configs()` entries and files next to the including file; `radiopyo/data/<x>` or a bare bundled name suggests `radiopyo/<x>`.
- [ ] A relative entry found only via the working directory logs a WARNING naming both the including file's folder and the CWD path used.
- [ ] Circular includes raise a `ConfigError` (still a `ValueError`) naming the cycle.
- [ ] `rp.bundled_configs()` returns `{"radiopyo/<posix path below data>": absolute Path}` for every `*.toml` under `PACKAGE_DATA_PATH`, recursively, sorted.
- [ ] Existing include tests pass (update the expected exception type where needed); new tests in `tests/parser/test_include_diagnostics.py`.

**Files to Touch:**
- `radiopyo/parser/exceptions.py`
- `radiopyo/parser/file_parser.py`
- `radiopyo/parser/utils.py`
- `radiopyo/parser/__init__.py`
- `radiopyo/__init__.py`
- `tests/parser/test_includes.py`
- `tests/parser/test_include_diagnostics.py`

**Verification Step:**

Run:
```bash
C:/Users/rtonneau/AppData/Local/miniconda3/envs/radiobio_env/python.exe -m pytest tests/parser -q && C:/Users/rtonneau/AppData/Local/miniconda3/envs/radiobio_env/python.exe -m ruff check radiopyo && C:/Users/rtonneau/AppData/Local/miniconda3/envs/radiobio_env/python.exe -m mypy radiopyo
```

Expected:
All parser tests pass; ruff "All checks passed!"; mypy "Success: no issues found".

**Notes:**

Anchoring currently rewrites the entry in `anchor_includes` and loses what was written: keep the original entry and the candidates tried (e.g. resolve in one place, `include_path` returning the path plus the tried list) so the error can report them. Keep `radiopyo/` the only bundled prefix (no silent alias). Use the radiobio_env python.exe path above.
