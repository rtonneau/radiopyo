# Ticket 02: per-file-parse-errors

**Model:** medium
**Effort:** medium

**Acceptance Criteria:**
- [ ] `TOMLFileParser.parse` catches `tomllib.TOMLDecodeError` and raises `ConfigSyntaxError` naming the file (and the include chain when it is an included file), original chained as `__cause__`; an invalid-escape message gets the hint: use single quotes ('C:\path\x.toml') or forward slashes.
- [ ] Every reaction string of a file (all reaction sections) is parsed with `RPARSER` when the file is read; an invalid one raises `ConfigError` naming the file, the section (`[[reactions.k_reaction]]`) and the string.
- [ ] Tests cover a syntax error in a main file and in an included file, a Windows-path escape, and a bad reaction string in an included file.

**Files to Touch:**
- `radiopyo/parser/toml.py`
- `radiopyo/parser/file_parser.py`
- `tests/parser/test_parse_errors.py`

**Verification Step:**

Run:
```bash
C:/Users/rtonneau/AppData/Local/miniconda3/envs/radiobio_env/python.exe -m pytest tests/parser -q && C:/Users/rtonneau/AppData/Local/miniconda3/envs/radiobio_env/python.exe -m ruff check radiopyo && C:/Users/rtonneau/AppData/Local/miniconda3/envs/radiobio_env/python.exe -m mypy radiopyo
```

Expected:
All parser tests pass; ruff and mypy report zero errors.

**Notes:**

The include chain is known in `resolve_includes`: wrap the parse call there (or pass the chain to the parser) so a nested error shows `main.toml -> a.toml`.
