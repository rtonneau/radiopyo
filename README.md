<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/assets/logo-dark.svg">
    <source media="(prefers-color-scheme: light)" srcset="docs/assets/logo.svg">
    <img src="docs/assets/logo.svg" alt="radiopyo" width="360">
  </picture>
</p>

<p align="center">
  <strong>Water radiolysis and radiobiological chemistry under conventional and ultra-high-dose-rate irradiation.</strong>
</p>

<p align="center">
  Dose rate and pulse timing change how long reactive species coexist, so the same dose
  can lead to different chemistry. radiopyo follows that chemistry in time.
</p>

<p align="center">
  <a href="https://github.com/rtonneau/radiopyo/actions/workflows/tests.yml"><img src="https://github.com/rtonneau/radiopyo/actions/workflows/tests.yml/badge.svg" alt="Tests"></a>
  <img src="https://img.shields.io/badge/python-%E2%89%A53.11-12665E" alt="Python 3.11 or later">
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-LGPL--3.0-12665E" alt="License: LGPL-3.0"></a>
  <a href="https://doi.org/10.1088/1361-6560/ae1fcc"><img src="https://img.shields.io/badge/DOI-10.1088%2F1361--6560%2Fae1fcc-BB583D" alt="DOI 10.1088/1361-6560/ae1fcc"></a>
</p>

## What you can study

radiopyo is a **homogeneous chemical-kinetics model**. Primary radiolytic yields (G-values)
feed a system of rate equations for the concentrations of the species, which
[scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)
integrates (LSODA by default, with an analytic Jacobian). It does not simulate particle
tracks or the spatial structure of the energy deposition. It started as a port of
R. Labarbe's MatLab model.

- **Dose rate and beam timing.** Constant, single-pulse, pulsed and macro-pulsed beams
  (trains of micro pulses), with a maximum dose or a number of pulses. Pulsed runs are
  cut exactly at the beam ON/OFF transitions.
- **Chemical interactions.** Water radiolysis chemistry plus a biological layer (oxygen,
  organic radicals R· and ROO·, catalase, Fenton chemistry). Rate-constant, Arrhenius,
  acid/base, equilibrium and Michaelis–Menten reactions, all written in TOML.
- **Oxygen and pH.** Initial concentrations (e.g. dissolved O₂) are inputs. The pH sets the
  initial H₃O⁺/OH⁻ and can be held constant (`force_cst_pH`).
- **LET.** A beam LET in 0.13–150 keV/µm sets the G-values of e⁻aq, OH·, H·, H₂, H₂O₂ and
  H₃O⁺ from the [Boscolo et al. (2020)](https://doi.org/10.3390/ijms21020424) fit. An
  optional dose-rate correction of these G-values exists, but its validity range is not
  checked.
- **Temperature.** Rate constants can depend on temperature (Arrhenius or polynomial).
- **Analysis.** Concentrations over time, per-reaction rates, production and consumption
  of each species, time-integrated concentrations, and HDF5 storage of runs.

## Quick start

radiopyo needs **Python 3.11 or later**; the dependencies (numpy, scipy, pandas, lark,
h5py, more-itertools) are installed automatically.

```bash
pip install git+https://github.com/rtonneau/radiopyo.git
```

Run the bundled Labarbe configuration (pure-water radiolysis plus the biological layer,
with a 10 Gy constant beam at 1e7 Gy/s) for 10 s:

```python
import radiopyo as rp

sim = rp.UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
res = sim.run(t_span=[1e-7, 10])   # [s]

df = res.to_pandas()               # concentrations [µmol/L], one column per species
print(df[["O2", "H2O2"]].iloc[-1]) # state at t = 10 s
```

`res` is an `ODEResult`. `to_pandas()` gives a DataFrame indexed by time [s], one column
per species (fixed species included), in µmol/L. The same object gives the
time-integrated concentrations (`res.integrate_species()`) and the contribution of each
reaction (`res.integrate_production_of("H2O2")`, `res.integrate_consumption_of("O2")`,
`res.reconstruct_reactions()`).

> [!NOTE]
> The radiolytic yield of the first 1e-7 s of a beam ON phase is added analytically to
> the initial concentrations, so the result starts at `t_span[0] + 1e-7`.

## Irradiation and chemistry

To compare dose rates, deliver the same dose with two different beams and look at the
chemistry afterwards. Here 10 Gy is delivered at 0.1 Gy/s (CONV) and in 5 pulses of
2 µs at 1e6 Gy/s (UHDR), and the oxygen consumed is read 200 s after the start:

```python
import radiopyo as rp

beams = {
    "CONV": rp.ConstantBeam(0.1, max_dose=10),          # 0.1 Gy/s for 100 s
    "UHDR": rp.PulsedBeam.from_dose_rate_peak(
        1e6,          # peak dose rate [Gy/s]
        0.1,          # period [s], i.e. 10 Hz
        2e-6,         # on_time [s]: 2 Gy per pulse
        max_dose=10,  # [Gy], by keyword (the 4th positional argument is n_pulse)
    ),
}

for name, beam in beams.items():
    sim = rp.UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
    sim.set_beam(beam, LET=0.201)  # [keV/µm], sets the G-values
    res = sim.prepare_chunked_run(t_span=[0, 200]).run()
    o2 = res.to_pandas()["O2"]
    print(f"{name}: O2 consumed = {o2.iloc[0] - o2.iloc[-1]:.2f} µmol/L")
```

`prepare_chunked_run` splits the run into beam ON and OFF chunks. Each chunk is integrated
with a constant dose rate and can have its own maximum step (`max_step_size_on`,
`max_step_size_off`). Which yield goes up or down with the dose rate depends on the
chemistry, the oxygen content and the beam; the model gives no general rule.

Beam timing affects how the run is solved:

- A configured beam ON time below 1e-6 s raises `ShortBeamOnError` when the beam is built.
- A run-time ON span shorter than that (a last pulse cut by `max_dose`, a pulse clipped by
  `t_span`) is not integrated. Its radiolytic yield is added analytically, a warning is
  logged, and the span is listed in `res.short_on_spans`.
- If scipy fails on a chunk, `SolverFailedError` is raised instead of returning a
  truncated result.
- Long pulsed runs store every solver step. `points_per_chunk_on` / `points_per_chunk_off`
  keep N points per chunk instead; integrals computed from such a thinned result are
  coarser and log a warning.
- `set_beam(beam)` without `LET` reuses the LET of the previously selected beam.

## Documentation

- **[Configuration reference](docs/configuration.md#configuring-simulation)**: includes,
  parameters, concentrations and every reaction type.
- **[Beam definitions](docs/configuration.md#possible-beam-definitions)**: constant, pulsed,
  single-pulse and macro-pulsed beams, LET and G-values.
- **[Running simulations](docs/configuration.md#ways-of-running-simulation)**: single and
  chunked runs, point thinning, solver failures, overshoot handling.
- **[Reaction model](docs/20230425_RPM_extended_doc.pdf)** (PDF): the equations and rate
  laws behind the ODE system.
- **[Bundled configurations](radiopyo/data/)**: the TOML files shipped with the package.
- **[Basic usage notebook](examples/notebook_basic_usage.ipynb)**: a worked example.
- **[Paper reproduction notebook](examples/tonneau_2025_pmb/tonneau_2025_pmb.ipynb)**:
  the main figures of Tonneau et al. 2025 (O2, H2O2 and pH time evolution, ΔO2 and net
  H2O2 yield maps over dose rate and pO2, gross production / consumption bars) with a
  comparison to the values quoted in the paper. It loads the recorded grid by default
  (`RECOMPUTE = True` re-runs it, about 10 min); needs `pip install "radiopyo[examples]"`.
- **[GF(DR) fit comparison notebook](examples/tonneau_2025_pmb/gf_fit_comparison.ipynb)**:
  the same grid run with the GF(DR) fit printed in the paper, compared with the shipped fit
  and with the values quoted in the paper.
- **[Fenton showcase](docs/issues/fenton/Fenton_Fixed_Showcase.ipynb)** (notebook): a Fenton
  runaway case and maps over dose per pulse and average dose rate.
- **[Changelog](CHANGELOG.md)**: changes between versions.

The bundled configurations are `advanced_water_radiolysis.toml` (pure water,
`rp.PURE_WATER_CONFIG_FILE`), `reactions_biology.toml` (the biological layer, to be used
with the former) and `config_Labarbe.toml` (both plus a constant beam,
`rp.LABARBE_CONFIG_FILE`). Your own file can include them and override single values.

<details>
<summary><strong>Configuration sections at a glance</strong></summary>

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
| `[[reactions.enzymatic]]` | Michaelis–Menten reaction (`k_value`, `k_micha`) |
| `[[beam.constant]]` | Constant dose rate, optional `max_dose` |
| `[[beam.pulsed]]` | Regular pulses (`dose_rate_peak`, `on_time`, `period`) |
| `[[beam.singlepulse]]` | One pulse of duration `on_time` |
| `[[beam.macro_pulsed]]` | Trains of micro pulses |

</details>

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

## Development

Clone the repository and install it as an editable package, with the test and type-check
tools:

```bash
git clone https://github.com/rtonneau/radiopyo.git && cd radiopyo
pip install -e ".[dev]" ruff
ruff check radiopyo && mypy radiopyo && python -m pytest
```

Ruff and mypy must report zero errors. Optional: `pip install "radiopyo[tqdm]"` for a
progress bar. To update a regular install, run
`pip install --upgrade git+https://github.com/rtonneau/radiopyo.git`; for an editable
install, `git pull`.

`tests/golden/` compares the outputs with recorded reference data, including the grid of
the paper reproduction (`tests/golden/data/paper_tonneau_2025.npz`). Re-record it only for
an approved result-changing fix:

```bash
PYTHONHASHSEED=0 python tests/golden/record_golden.py --force --only paper --workers 8
```

## Acknowledgement

This module was originally inspired by R. Labarbe's
[validated](https://doi.org/10.1016/j.radonc.2020.06.001) MatLab code.

## License

radiopyo is released under the [LGPL-3.0](LICENSE) license.
