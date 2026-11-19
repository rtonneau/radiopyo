# Updating radiopyo on a Windows desktop
<!-- radiopyo-version: 1.1.0 (scripts/bump_doc_version.py rewrites this version below at each release) -->

This guide updates an existing radiopyo (0.x, often installed from the former GitLab
repository `gitlab.unamur.be`) to the current release, **v1.1.0**, on Windows with
Anaconda or Miniconda. Every command is typed in the **Anaconda Prompt** (Start menu →
"Anaconda Prompt").

radiopyo now lives on GitHub: <https://github.com/rtonneau/radiopyo> (public, no login
needed). The GitHub history continues the GitLab one.

What changed between your version and v1.1.0 is listed in the [changelog](../CHANGELOG.md).

## 1. Find out how the old version was installed

Activate the environment that holds the old radiopyo (replace `radiobio` with its name;
`conda env list` lists them), then look at how pip sees it:

```bat
conda activate radiobio
python --version
pip freeze | findstr /i radiopyo
```

Note the Python version, then compare the `pip freeze` line with this table (ignore a
`radiopyo-rs` / `radiopyo_rs` line: that is the optional Rust extension of step 6):

| `pip freeze` line | Case | Go to |
|---|---|---|
| `radiopyo @ git+https://gitlab.unamur.be/...` (or `github.com`) | installed with pip from a git URL | [step 3a](#3a-installed-with-pip-from-a-git-url) |
| `-e git+...#egg=radiopyo`, or `-e C:\...` (after a line `# Editable install ...`) | editable install of a local clone | [step 3b](#3b-editable-install-of-a-local-clone) |
| `radiopyo==0.x.y` or `radiopyo @ file:///...whl` | installed from a wheel or an archive | [step 3a](#3a-installed-with-pip-from-a-git-url) |
| nothing | radiopyo is not in this environment | check another one with `conda env list` |

## 2. Choose the environment

radiopyo v1.1.0 needs **Python 3.11 or later**. Old versions accepted Python 3.9.

- **`python --version` shows 3.11 or later:** you can update in place, in the same
  environment.
- **It shows 3.9 or 3.10:** create a new environment. Keep the old one if you want to
  re-run old scripts with the old version.

  ```bat
  conda create -n radiopyo python=3.11
  conda activate radiopyo
  ```

  In the new environment, skip the uninstall of step 3a; for an editable clone (step 3b)
  run everything but the uninstall.

pip installs from a git URL, so git must be available:

```bat
git --version
```

If Windows answers `'git' is not recognized ...`, install it in the active environment
with `conda install git` (or install [Git for Windows](https://git-scm.com/download/win)
and open a new Anaconda Prompt).

## 3a. Installed with pip from a git URL

Remove the old version, then install v1.1.0 from GitHub:

```bat
pip uninstall -y radiopyo
pip install "radiopyo[tqdm,examples] @ git+https://github.com/rtonneau/radiopyo.git@v1.1.0"
```

The extras are optional: `tqdm` gives the progress bar (it used to be installed
automatically), `examples` installs matplotlib for the notebooks. `pip install
"radiopyo @ git+https://github.com/rtonneau/radiopyo.git@v1.1.0"` installs radiopyo alone.

To follow the development version instead of the release, drop `@v1.1.0` from the URL.

## 3b. Editable install of a local clone

The `-e C:\...` path of step 1 is the clone folder. Go there and check for local changes
first:

```bat
cd /d C:\path\to\radiopyo
git status
```

If `git status` lists modified files you want to keep, commit them (`git commit -am "my
changes"`) or set them aside (`git stash`). Then point the clone at GitHub and fetch the
release:

```bat
git remote -v
git remote set-url origin https://github.com/rtonneau/radiopyo.git
git fetch origin --tags
git checkout v1.1.0
```

`git remote -v` shows where `origin` pointed before (the GitLab URL for an old clone).
`git checkout v1.1.0` puts the clone exactly on the release (git reports a "detached
HEAD": expected). To follow the development version instead, run `git checkout main` then
`git pull origin main`.

Reinstall the package so pip picks up the new dependencies:

```bat
pip uninstall -y radiopyo
pip install -e ".[tqdm,examples]"
```

## 4. Check the update

In the same Anaconda Prompt:

```bat
python -c "import importlib.metadata as m, radiopyo; print(m.version('radiopyo'), radiopyo.__file__)"
```

It must print `1.1.0` and a path inside the environment (`...\envs\<name>\Lib\site-packages\radiopyo\...`),
or inside your clone for an editable install. Then run the Quick start example: save these
lines as `check_radiopyo.py` in any folder **other than** a radiopyo clone and run
`python check_radiopyo.py`.

```python
import radiopyo as rp

sim = rp.UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
res = sim.run(t_span=[1e-7, 10])   # [s]
print(res.to_pandas()[["O2", "H2O2"]].iloc[-1])
```

It prints the O2 and H2O2 concentrations (µmol/L) at t = 10 s, after a few INFO log lines.

## 5. Jupyter notebooks

A notebook runs the radiopyo of its **kernel**, not of the Anaconda Prompt. If you created
a new environment in step 2, register it as a kernel:

```bat
pip install ipykernel
python -m ipykernel install --user --name radiopyo --display-name "Python (radiopyo)"
```

Then select "Python (radiopyo)" in Jupyter or VS Code (Kernel → Change kernel). In a cell,
`import sys; print(sys.executable)` shows which environment the kernel uses.
`jupyter kernelspec list` lists the registered kernels; `jupyter kernelspec uninstall
<name>` removes an old one.

## 6. Rust backend (optional)

v1.1.0 can solve a simulation with radkin, a Rust port of radiopyo, through the
`radiopyo-rs` extension (LSODA only). On Windows, install the
wheel of the release (64-bit Python 3.11 or later, no Rust toolchain needed):

```bat
pip install https://github.com/rtonneau/radiopyo/releases/download/v1.1.0/radiopyo_rs-1.1.0-cp311-abi3-win_amd64.whl
python -c "import radiopyo as rp; print(rp.rust.is_available())"
```

It prints `True` when the extension loads. Then pass `backend="rust"`:

```python
res = sim.run(t_span=(1e-7, 10), backend="rust")
res = sim.prepare_chunked_run(t_span=(0, 1e-2)).run(backend="rust")
```

See the [Rust backend notebook](../examples/notebook_rust_backend.ipynb) for a worked
example. The extension's version must match radiopyo's: when you update radiopyo later,
install the wheel of the same release.

## Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `'git' is not recognized as an internal or external command` | git is missing: `conda install git` in the environment (step 2). |
| `ERROR: Package 'radiopyo' requires a different Python: 3.9.x not in '>=3.11'` | The environment is too old: create a Python 3.11 environment (step 2). |
| Step 4 still prints an old version | Several copies are installed. Run `pip uninstall radiopyo` until pip says it is not installed, then install again. `where python` shows which Python the prompt uses; check that the right environment is active. |
| Step 4 prints a path outside the environment | Python imports a `radiopyo` folder from the current directory (an old clone): run from another folder. |
| A notebook still uses the old version | Its kernel points at the old environment (step 5). Restart the kernel after any install. |
| `ModuleNotFoundError: No module named 'tqdm'` | `tqdm` is now an extra: `pip install tqdm`. |
| `ImportError: The Rust backend needs the optional 'radiopyo-rs' extension` | The wheel is not installed in this environment (step 6). |
| `radiopyo_rs-...whl is not a supported wheel on this platform` | Python is 32-bit or older than 3.11: use a 64-bit Python 3.11 environment. |
| An old script or config fails (`.ron` config, `force_positive`, Arrhenius `Ea`, short beam ON time) | Behaviour changed between versions: see the [changelog](../CHANGELOG.md), sections "Removed" and "Changed". HDF5 result files written by old versions still load. |
