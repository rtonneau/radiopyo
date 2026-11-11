# radiopyo: Python Module dedicated to Radiobiology

![version](https://img.shields.io/badge/version-0.6.0-blue)
![License](https://img.shields.io/badge/License-LGPL3-radiopyo.svg)
![maintained](https://img.shields.io/badge/Maintained%3F-yes-green.svg)
![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)

This python3 module allows you to perform radiobiology related simulations. The embedded
science (all Physics, Chemistry & Biology) has the purpose to solve systems of Ordinary
Differential Equations (a.k.a. ODE) based on equations provided via an input file.
The default solver is the LSODA solver from the good
[scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)
module. The solver uses the analytic Jacobian by default (`use_jac=True`, honoured by
LSODA, BDF and Radau; the other methods ignore it). Pass `use_jac=False` to `run` to
fall back on the finite-difference Jacobian of scipy.

## Remarks for developer

If you intend to join this dev. adventure, first of all, thank you for your time and
support. Second, please be aware that this code base uses the fast
[RUFF](https://github.com/astral-sh/ruff) linter and comply with the
[Mypy](https://mypy.readthedocs.io/en/stable/index.html) static type checker.

ps: Needless to say that Python2 is forbidden...

Tests live in `tests/` (`python -m pytest`), including "golden" outputs recorded from
the code before the refactoring. Developers install the extra tools with
`pip install -e ".[dev]"` (pytest, mypy, pandas-stubs), then run
`ruff check radiopyo` and `mypy radiopyo`.

## Getting Started

This documentation assumes you already have or know how to install a new python3
environment on your local machine. If it is not the case, we suggest to check
[this](https://www.anaconda.com/) out.

Python >= 3.11 is required. The dependencies of this module (installed automatically)
are:

- scipy
- numpy
- pandas
- lark
- more-itertools
- h5py (storage of results in HDF5 files)

For a nice progress bar, install the optional `tqdm` extra
(`pip install "radiopyo[tqdm]"`). JupyterLab users should also install
`jupyterlab-widgets` and `ipywidgets`.

## Installation

1. Via pip + git
   1. Run:

      ```bash
      pip install git+https://github.com/rtonneau/radiopyo.git
      ```

   2. You are good to go.

2. Install as an editable package for development
   1. Clone this repo
   2. In the package folder, run the following:

      ```bash
      pip install -e ./          # or pip install -e ".[dev]" to also get the dev tools
      ```

   3. You are good to go.

## Usage

   ```python
   import radiopyo as rp

   # Simulation handler
   # rp.LABARBE_CONFIG_FILE contains environment and beam definition
   sim = rp.UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
   res = sim.run(t_span=[1e-7, 10]) #Run simulation from 1e-7 s to 10 s.

   # Extract results to a nice pandas DataFrame.
   df = res.to_pandas()

   # Time integrate species concentrations.
   s = res.integrate_species() # pandas.Series
   ```

## Configuring Simulation

For the moment, the only ways to load simulation configurations is via
[toml](https://toml.io/en/) file.
The TOML configuration file consists of several sections:

### \[includes\]

This section enable the user to reference other config file to read.
All included files will have a lower precedence i.e. will be erased by new definitions.
It is then, for example, possible to modify a rate constant of a given reaction defined
in an included file.

```toml
[includes]
toml = [
    "radiopyo/advanced_water_radiolysis.toml",
    "radiopyo/reactions_biology.toml",
]
```

**Where the files are looked for**

- An entry starting with exactly "radiopyo/" is referring to internal data
  (radiopyo/radiopyo/data) i.e. provided in the package:
  "radiopyo/reactions_biology.toml" is the "reactions_biology.toml" shipped with the
  package. Any other entry is a plain path, even if it starts with "radiopyo"
  (e.g. "radiopyo_extra.toml").
- A relative path is resolved against the folder of the file containing the entry
  (so nested includes keep working wherever the program is started from). If the file
  is not there, the current working directory is tried.
- An absolute path towards a valid configuration file is always possible.

**Precedence**

From the lowest to the highest: the included files, in the order of the list (a file
listed later overrides the ones listed before it), then the file including them. Each
included file may itself include other files.

**How the sections are merged**

- `[parameters]`, `[concentrations.fixed]`, `[concentrations.initial]` and
  `[reactions.radiolytic]` are merged entry by entry: a name defined again takes the new
  value, the others are kept.
- Reactions (`k_reaction`, `constant`, `arrhenius`, `acid_base`, `equilibrium` and
  `enzymatic`) are identified by their species on each side of the arrow, whatever their
  order: "B + A -> C" redefines "A + B -> C", but "C -> A + B" is another reaction.
  A redefinition updates the reaction field by field: only what is written is replaced,
  the rest is kept. For example, writing `k_forward` alone for an acid/base reaction keeps
  its `k_backward` and `pKa`, and writing `k_value` alone for an enzymatic reaction keeps
  its `k_micha`. The only exception is an `equilibrium` giving at least two of
  `k_forward`, `k_backward` and `k_equilibrium`: it is then fully defined by the new
  entry and the constants of the included one are dropped. The values are never
  rescaled (for "2 A -> B" the value to write is still 2k). Redefined reactions keep their
  place, new ones are appended.
- Beams: see "Possible beam definitions" below.

Files defined at the package level are the following:

- advanced_water_radiolysis.toml
   Contains all reactions related to the radiolysis of pure water.
- reactions_biology.toml
   Adds the biological part (fixed species, initial O2, radiolytic yield of R_r and the
   related reactions). It has no `[includes]` section: use it together with
   'advanced_water_radiolysis.toml'.
- config_Labarbe.toml
   Includes both 'advanced_water_radiolysis.toml' and 'reactions_biology.toml' and add a
   constant beam matching the one used in the excellent paper of R. Labarbe.

In Python, `rp.PURE_WATER_CONFIG_FILE` and `rp.LABARBE_CONFIG_FILE` give the path of
'advanced_water_radiolysis.toml' and 'config_Labarbe.toml'.

### \[parameters\]

Contains the main parameters. All are optional:

| Key | Default | Meaning |
| --- | --- | --- |
| `pH` | 7 | Initial pH |
| `force_cst_pH` | false | Keep the pH constant during the simulation |
| `temperature` | 25 | Simulation temperature [degC] |
| `variable_GValues` | false | Apply a dose-rate dependent correction to the LET dependent G-values (see "Possible beam definitions") |
| `species_order` | (none) | Explicit order of the species in the state vector |

```toml
[parameters]
pH = 7.4
temperature = 37
```

The temperature is read by `rp.UnitCell` and can be changed afterwards with
`sim.set_temperature(37)` (or `sim.T = 37`).

### \[concentrations.fixed\]

Contains the concentrations of species supposed to be constant.
=> If a species is recorded in this section, its concentration will be forced constant.
For example:

```toml
[concentrations.fixed] # Unit is [mol]/[l]
catalase = 8.0e-8
```

### \[concentrations.initial\]

Contains the initial concentrations of the species.
=> If a species is not recorded in this section, its initial concentration will be
assumed as 0.0 mol/l.
For example:

```toml
[concentrations.initial] # Unit is [mol]/[l]
O2 = 50e-6
```

### \[reactions.radiolytic\]

Contains the G-values of the species created by the radiation
source. The unit is [radical/100eV/incident particle].
For example:

```toml
[reactions.radiolytic]
e_aq = 2.8
OH_r = 2.8
```

### \[\[reactions.acid_base\]\]

Contains all the acid/base reactions involved in the model. Each
entry must define 4 fields: "reaction" (the acid on the left of the "<->" arrow, its
conjugated base and "H3O_plus" on the right, or "OH_minus" as a reactant), "pKa",
"k_forward" and "k_backward". A "k_value" is rejected: it is ambiguous for a
forward/backward pair.
For example:

```toml
[[reactions.acid_base]]
reaction = "H2O2 <-> H3O_plus + HO2_minus"
pKa = 11.65
k_forward = 9.43e-2
k_backward = 5.02e10

```

ps: note the double "\[\[" and "\]\]", mandatory for each entry!

### \[\[reactions.equilibrium\]\]

Same as "acid_base" for a pair of forward/backward reactions that is not an acid/base
couple. The reaction string uses the "<->" arrow and 2 out of the 3 constants
"k_forward", "k_backward" and "k_equilibrium" must be given (the third one is deduced).

```toml
[[reactions.equilibrium]]
reaction = "H_r <-> e_aq + H3O_plus"
k_forward = 5.83
k_backward = 2.09e10
```

### \[\[reactions.k_reaction\]\]

Contains all the simple chemical reactions involved in the model.
Each list entry must defined 2 fields: "reaction" & "k_value".
For example:

```toml
[[reactions.k_reaction]]  # O_r_minus + O2_r_minus -> 2 OH_minus + O2
reaction = "O_r_minus + O2_r_minus -> 2 OH_minus + O2"
k_value = 6.0e8
```

Charges must be specified either by '_minus' or '_plus'. For multiple charge, just
multiply the specifier e.g. $Fe^{2+}$ -> Fe_plus_plus.
Reactants and products must only be separated by the '->' string.

**Convention for n A -> B.** For a reaction with a stoichiometric coefficient n on a
reactant (e.g. "2 OH_r -> H2O2"), the "k_value" to write is **n times the rate constant
k** (the value is divided by n when the reaction is built). This applies to a single
number only: a temperature polynomial, an Arrhenius `A` or an equilibrium constant is
used as given.

Blanks/spaces for separation are allowed but not mandatory.

ps: note the double "\[\[" and "\]\]", mandatory for each entry!

### \[\[reactions.arrhenius\]\]

Same as "k_reaction" for the reaction string, but the rate constant depends on the
temperature of the simulation (`temperature` in `[parameters]`, in degC). Each list entry
must define 3 fields: "reaction", "A" and "Ea". The constant is

$k(T) = A \exp\left(-\frac{E_a}{R\,T}\right)$

with `A` the pre-exponential factor (same unit as a "k_value", e.g. L/mol/s for a
bimolecular reaction), `Ea` the activation energy in **kJ/mol**, $R$ = 8.314 J/mol/K and
$T$ the temperature in K. For example (k = 3.13e6 L/mol/s at 25 degC):

```toml
[[reactions.arrhenius]]
reaction = "H_r + OH_r -> H2O"
A = 1e10    # [L/mol/s]
Ea = 20     # [kJ/mol]
```

Before this was documented, `Ea` was read as $E_a/R$ in K, i.e. $k = A \exp(-E_a/T)$:
multiply an old value by 8.314e-3 to get kJ/mol.

ps: note the double "\[\[" and "\]\]", mandatory for each entry!

### \[\[reactions.enzymatic\]\]

Contains all the Enzymatic reactions involved in the model. The
enzymatic reactions are modelled with the Michaelis-Menten kinetic. Each list
entry must defined 3 fields: "reaction", "k_value" & "k_micha".
For example:

```toml
[[reactions.enzymatic]]
reaction =" 2H2O2 -- catalase >> O2 + 2 H2O"
k_value = 6.62e7 #  [1/s]
k_micha = 1.1 #  [mol/l]

```

The same grammar as for "k_reactions" applies. For enzymatic reactions, the enzyme name
must be specified in the operator arrow (-- 'enzyme name' >>). It is a plain name, in any
case (`catalase`, `SOD`, `GPx`), never read as a chemical formula. As for k-reactions,
the rate of an enzymatic reaction is multiplied, for each of its "negative" products, by
c / (c + 1e-12 mol/L) (c = concentration of the product): it is halved at 1e-12 mol/L and
vanishes as the product is depleted. Rate laws always read the clipped concentration
max(c, 0), so a species that overshoots below 0 is never consumed.

ps: note the double "\[\[" and "\]\]", mandatory for each entry!

### Possible beam definitions

All of the following are valid beam definitions:

```toml
[[beam.constant]]
label          = "constant beam example"
dose_rate      = 0.25  # Average dose rate [Gy/s]
max_dose       = 7.5   # Max dose to deliver [Gy], after that, beam is OFF.
LET            = 0.201 # [keV/µm]
```

```toml
[[beam.pulsed]]
label          = "pulse beam example 1"
dose_rate_peak = 319   # Dose rate during beam ON [Gy/s]
period         = 121   # Beam period repetition [s]
on_time        = 0.1   # Time duration for beam ON [s]
n_pulse        = 8     # Number of pulses
start_time     = 30    # Time of first pulse [s]
LET            = 0.201 # [keV/µm]
```

```toml
[[beam.pulsed]]
label          = "pulse beam example 2"
dose_rate_peak = 125   # Dose rate during beam ON [Gy/s]
on_time        = 8e-3  # Time duration for beam ON [s]
period         = 20e-3 # Beam period repetition [s]
max_dose       = 1     # Max dose to deliver [Gy]
LET            = 7.5   # [keV/µm]
```

```toml
[[beam.pulsed]]
label          = "pulse beam example 3"
dose_rate_peak = 5.7e5 # Peak dose rate for pulsed beam [Gy/s]
on_time        = 2e-6  # Time [seconds] (used only for Pulsed beam)
frequency      = 90    # Beam frequency [Hz]
max_dose       = 10    # Max dose to deliver [Gy] (can also be 'n_pulse')
LET            = 3.2   # [keV/µm]
```

```toml
[[beam.singlepulse]]
label          = "singlepulse example"
dose_rate      = 2.2e6 # Average dose rate (for singlepulse beam, same as peak) [Gy/s]
on_time        = 3e-6  # Time duration for beam ON [s]
LET            = 0.201 # [keV/µm]
```

```toml
[[beam.macro_pulsed]]   # trains of micro pulses
label            = "macro pulsed example"
dose_rate_peak   = 1e6   # Dose rate during a micro pulse [Gy/s]
period           = 0.2   # Macro period (time between two train starts) [s]
on_time          = 0.01  # Duration of a train [s]
micro_period     = 1e-3  # Period of the micro pulses inside a train [s]
micro_on_time    = 2e-6  # Duration of a micro pulse [s]
max_dose         = 10    # Max dose to deliver [Gy], optional
start_time       = 0.05  # Time of the first train [s], optional
micro_start_time = 0     # Start of the first micro pulse inside a train [s], optional
LET              = 3.2   # [keV/µm]
```

`frequency` and `micro_frequency` [Hz] can replace `period` and `micro_period` (also for
`pulsed`). A `pulsed` beam can also be defined with `pulse_dose` [Gy] (dose per pulse)
instead of two of `dose_rate_peak`, `on_time`, `pulse_dose`.

ps: note the double "\[\[" and "\]\]", mandatory for each entry!

The 'label' keyword is not mandatory. If not specify, the 'default' label will be given.
However, be aware that any attempt to create another unlabelled (or labelled 'default')
beam will erase previous definition. Additionally, be aware that the "default" beam is
never included/imported from files specified in the [includes] section.

LET values will redefine the G-values of the following species:

- e_aq
- OH_r
- H_r
- H2
- H2O2
- H3O_plus (see below)

G-values are interpolated according to the work of [D.
Boscolo](https://doi.org/10.3390/ijms21020424)
Therefore, LET values must be in the range: [0.13, 150] keV/µm

Note on H3O_plus: its G-value is taken **equal to the one of e_aq** (charge balance of
the primary yields). The `H3O_plus` column of the Boscolo data file is not used.

With `variable_GValues = true` in `[parameters]`, the G-values are additionally scaled by
`0.97207704 * dose_rate_peak ** -0.10972` (fit of the dose-rate dependence, with the peak
dose rate of the beam in Gy/s). The code does not check the dose-rate range of this fit
(and the factor is infinite at 0 Gy/s): stay within the range of the data it was fitted on.

## Ways of Running Simulation

The module currently supports 2 ways of running simulations:

- run in one single block

   ```python
   import radiopyo as rp

   # Simulation handler
   sim = rp.UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
   res = sim.run(t_span=[1e-7, 10]) #Run simulation.

   ```

   The initial radiolytic yield of the first 1e-7 s of a beam ON phase is added
   analytically to the initial concentrations: the result starts at
   `t_span[0] + 1e-7`.

- run in several chunks (useful for pulsed beams)

   ```python
   import radiopyo as rp

   # Beam Definition
   dose_rate_peak  = 1000 # Gy/s
   max_dose        = 10   # Gy
   period          = 1e-3 # seconds
   on_time         = 5e-6 # seconds
   beam = rp.PulsedBeam.from_dose_rate_peak(dose_rate_peak, period, on_time, max_dose=max_dose)

   # Simulation handler
   sim = rp.UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
   sim.set_beam(beam)

   # Better to adapt simulation settings during ON and OFF time.
   res = sim.prepare_chunked_run(t_span=[1e-7, 10],
                                 max_step_size_on=1e-8,  #seconds
                                 max_step_size_off=0.01, #seconds
                                 ).run()

   ```

   The chunks are cut exactly where the beam switches (pulse start and end, time of
   the maximum dose) and the solver holds the dose rate of each chunk constant. A
   pulse must be longer than `UnitCell.T_ZERO` (1e-7 s), otherwise a `ValueError`
   ("Too Short beam pulse") is raised.

   If `scipy.integrate.solve_ivp` fails on a chunk, `SolverFailedError`
   (`radiopyo.simulation.exceptions`) is raised with the span, the last time reached and
   the message of scipy, instead of returning a truncated solution.

   A species slightly below 0 (a solver overshoot) is read as 0 by every
   reaction; only its own balance is continued linearly from 0 (overshoot relaxation,
   `.work/adr/0007-overshoot-relaxation.md`), which keeps LSODA from stalling on stiff
   radicals such as OH_r in a high-RH medium.

## Acknowledgement

- This module was originally inspired by R. Labarbe
  [validated](https://doi.org/10.1016/j.radonc.2020.06.001) MatLab code.
