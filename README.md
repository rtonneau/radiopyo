<p align="center"><img src="docs/logo.svg" alt="radiopyo" width="420"></p>

[![Tests](https://github.com/rtonneau/radiopyo/actions/workflows/tests.yml/badge.svg)](https://github.com/rtonneau/radiopyo/actions/workflows/tests.yml)
![version](https://img.shields.io/badge/version-0.6.0-blue)
![Python](https://img.shields.io/badge/Python-%3E%3D3.11-blue)
![License](https://img.shields.io/badge/License-LGPL--3.0-blue)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
![mypy](https://img.shields.io/badge/mypy-checked-blue)
[![DOI](https://img.shields.io/badge/DOI-10.1088%2F1361--6560%2Fae1fcc-blue)](https://doi.org/10.1088/1361-6560/ae1fcc)

radiopyo simulates water radiolysis and radiobiology chemistry under UHDR (FLASH) and
conventional irradiation. It builds a system of ordinary differential equations from a
[TOML](https://toml.io/en/) configuration file and solves it with
[scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)
(LSODA by default, with an analytic Jacobian).

## Installation

Python >= 3.11 is required; dependencies are installed automatically.

```bash
pip install git+https://github.com/rtonneau/radiopyo.git
```

For development, clone the repository and install it as an editable package:

```bash
git clone https://github.com/rtonneau/radiopyo.git && cd radiopyo && pip install -e ./
```

Optional: `pip install "radiopyo[tqdm]"` for a progress bar.

## Update

```bash
pip install --upgrade git+https://github.com/rtonneau/radiopyo.git
```

For an editable install, run `git pull` instead.

## Usage

Single run:

```python
import radiopyo as rp

sim = rp.UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
res = sim.run(t_span=[1e-7, 10])  # from 1e-7 s to 10 s

df = res.to_pandas()          # results as a pandas DataFrame
s = res.integrate_species()   # time-integrated concentrations (pandas.Series)
```

The radiolytic yield of the first 1e-7 s of a beam ON phase is added analytically to the
initial concentrations, so the result starts at `t_span[0] + 1e-7`. Pass `use_jac=False`
to `run` to replace the analytic Jacobian by scipy's finite differences.

Chunked run, useful for pulsed beams (separate step sizes for beam ON and OFF):

```python
import radiopyo as rp

beam = rp.PulsedBeam.from_dose_rate_peak(
    1000,         # dose_rate_peak [Gy/s]
    1e-3,         # period [s]
    5e-6,         # on_time [s]
    max_dose=10,  # [Gy], pass it by keyword
)

sim = rp.UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
sim.set_beam(beam)
res = sim.prepare_chunked_run(
    t_span=[1e-7, 10],
    max_step_size_on=1e-8,   # [s]
    max_step_size_off=0.01,  # [s]
).run()
```

A beam ON time shorter than 1e-6 s raises `ShortBeamOnError`. If scipy fails on a chunk,
`SolverFailedError` is raised instead of returning a truncated result.

## Configuration

A simulation is described by a TOML file. Full reference: [docs/configuration.md](docs/configuration.md).

| Section | Purpose |
| --- | --- |
| `[includes]` | Other config files to load first; the including file overrides them |
| `[parameters]` | `pH`, `force_cst_pH`, `temperature`, `variable_GValues`, `species_order` |
| `[concentrations.fixed]` | Species held at a constant concentration [mol/L] |
| `[concentrations.initial]` | Initial concentrations [mol/L] (default 0) |
| `[reactions.radiolytic]` | G-values of the species created by the radiation |
| `[[reactions.k_reaction]]` | Reaction with a rate constant `k_value` |
| `[[reactions.arrhenius]]` | Reaction with a temperature-dependent rate (`A`, `Ea` in kJ/mol) |
| `[[reactions.acid_base]]` | Acid/base pair with `pKa`, `k_forward`, `k_backward` |
| `[[reactions.equilibrium]]` | Forward/backward pair, 2 of 3 constants given |
| `[[reactions.enzymatic]]` | Michaelis-Menten reaction (`k_value`, `k_micha`) |

| Beam | Description |
| --- | --- |
| `[[beam.constant]]` | Constant dose rate, optional `max_dose` |
| `[[beam.pulsed]]` | Regular pulses (`dose_rate_peak`, `on_time`, `period`) |
| `[[beam.singlepulse]]` | One pulse of duration `on_time` |
| `[[beam.macro_pulsed]]` | Trains of micro pulses |

A beam `LET` [keV/um, 0.13 to 150] sets the G-values ([Boscolo et al.](https://doi.org/10.3390/ijms21020424)).

Bundled configs live in [radiopyo/data/](radiopyo/data/): `advanced_water_radiolysis.toml`
(pure water, `rp.PURE_WATER_CONFIG_FILE`), `reactions_biology.toml` (biological part, use
it with the former) and `config_Labarbe.toml` (both plus a constant beam,
`rp.LABARBE_CONFIG_FILE`).

## Development

```bash
pip install -e ".[dev]" ruff
ruff check radiopyo && mypy radiopyo && python -m pytest
```

Ruff and mypy must report zero errors. radiopyo is released under the LGPL-3.0 license.

## Citation

If you use radiopyo in your work, please cite:

R. Tonneau, R. Labarbe, S. Lucas and A.-C. Heuskin, "Decoding the chemistry of the FLASH
effect: a physicochemical model of dose-rate, pH and oxygen-dependent H₂O₂ production",
Phys. Med. Biol. 70 235021 (2025),
[link](https://iopscience.iop.org/article/10.1088/1361-6560/ae1fcc/meta),
DOI [10.1088/1361-6560/ae1fcc](https://doi.org/10.1088/1361-6560/ae1fcc).

```bibtex
@article{Tonneau_2025,
  title   = {Decoding the chemistry of the {FLASH} effect: a physicochemical model of dose-rate, {pH} and oxygen-dependent {H$_2$O$_2$} production},
  author  = {Tonneau, Romain and Labarbe, Rudi and Lucas, St{\'e}phane and Heuskin, Anne-Catherine},
  journal = {Physics in Medicine \& Biology},
  volume  = {70},
  number  = {23},
  pages   = {235021},
  year    = {2025},
  doi     = {10.1088/1361-6560/ae1fcc}
}
```

## Acknowledgement

- This module was originally inspired by R. Labarbe
  [validated](https://doi.org/10.1016/j.radonc.2020.06.001) MatLab code.
