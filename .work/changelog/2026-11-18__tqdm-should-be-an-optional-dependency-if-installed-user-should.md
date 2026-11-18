---
bump: minor
floor: minor
---
### Added

- The tqdm progress bar of chunked runs can be turned off: `SimChunksIter.run(progress=False)` for one run,
  `radiopyo.set_progress_bar(False)` for every later run (`radiopyo.progress_bar_enabled()` reads the setting).
  `progress=None` (default) follows the global setting, which is on by default; an explicit `progress=` wins.
  tqdm stays optional (the `tqdm` extra): without it no bar is shown and nothing fails.
