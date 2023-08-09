# radiopyo: Python Module dedicated to Radiobiology

![version](https://img.shields.io/badge/version-0.3.0-blue)
![License](https://img.shields.io/badge/License-LGPL3-radiopyo.svg)
![maintained](https://img.shields.io/badge/Maintained%3F-yes-green.svg)
![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v1.json)

This python3 module allows you to perform radiobiology related simulations. The embedded
science (all Physics, Chemistry & Biology) has the purpose to solve systems of Ordinary
Differential Equations (a.k.a. ODE) based on equations provided via an input file.
The default solver is the LSODA solver from the good
[scipy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html)
module.

## Remarks for developer

If you intend to join this dev. adventure, first of all, thank you for your time and
support. Second, please be aware that this code base uses the fast
[RUFF](https://github.com/astral-sh/ruff) linter and comply with the
[Mypy](https://mypy.readthedocs.io/en/stable/index.html) static type checker.

ps: Needless to say that Python2 is forbidden...  
pps: There is currently no code coverage (tests)... a big problem to solve.

## Getting Started

This documentation assumes you already have or know how to install a new python3
environment on your local machine. If it is not yest the suggest, I suggest to check
[this](https://www.anaconda.com/) out.

Currently the main dependencies of this modules are:

- scipy
- numpy
- pandas
- pandas-stubs
- lark
- more-itertools

For JupyterLab user, think to also install (to have nice progressbar):

- tqdm
- jupyterlab-widgets
- ipywidgets

## Installation

1. Via pip + git
   1. Run:

      ```bash
      pip install git+https://gitlab.unamur.be/rtonneau/radiopyo.git
      ```

   2. You are good to go.

2. Install as an editable package for development
   1. Clone this repo
   2. In the package folder, run the following:

      ```bash
      pip install -e ./
      ```

   3. You are good to go.

## Usage

   ```python
   import radiopyo as rp 

   # Simulation handler
   # rp.LABARBE_CONFIG_FILE contains environment and beam definitions
   sim = rp.UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
   res = sim.run(t_span=[1e-9, 10]) #Run simulation from 1e-9 s to 10 s.

   # Extract results to a nice pandas DataFrame.
   df = res.to_pandas() 

   # ---------- Same but change beam definition ----------
   peak_dose_rate  = 1000 # Gy/s (dose rate during ON time)
   max_dose = 10 # Gy
   period = 10e-3 # seconds
   on_time = 1e-3 # seconds (=> Duty cycle of 10%)
   beam = rp.PulsedBeam.from_peak_dose_rate(peak_dose_rate, period, on_time, max_dose)
   
   sim = rp.UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
   sim.set_beam(beam)
   
   # Better to adapt simulation settings during ON and OFF time.
   res = sim.prepare_chunked_run(t_span=[1e-9, 10],
                                 atol=1e-4, #Absolute tolerance of the ODE Solver
                                 rtol=1e-6, #Relative tolerance of the ODE Solver
                                 max_step_size_on=1e-8,  #seconds
                                 max_step_size_off=0.01, #seconds
                                 ).run()
   # /!\ For atol and rtol, the base unit for the ODE solver is [µmol/l]

   df = res.to_pandas() # pandas.DataFrame
   
   # Time integrate species concentrations.
   s = res.integrate_species() # pandas.Series
   ```

## Configuring Simulation

For the moment, the only ways to load simulation configurations is via
[ron](https://github.com/ron-rs/ron) or [toml](https://toml.io/en/) files.
TOML files are recommended.
The TOML configuration file consists of several sections:

### [includes]

This section enable the user to reference other config file to read.
All included files will have a lower precedence i.e. will be erased by new definitions.
It is then, for example, possible to modify a rate constant of a given reaction defined
in an included file.

```toml
toml = [
    "radiopyo/basic_water_radiolysis.toml",
    "radiopyo/reactions_biology.toml",
]
```

"radiopyo/" is referring to internal data (radiopyo/radiopyo/data) i.e. provided in the
package. Of course, each entry can be an absolute path towards a valid configuration
file. Files defined at the package level are the following:

- basic_water_radiolysis.toml  
   Contains all reactions related to the radiolysis of pure water.
- reactions_biology.toml  
   Includes 'basic_water_radiolysis.toml' and add biological reactions.
- config_Labarbe.toml  
   Includes both 'basic_water_radiolysis.toml' and 'reactions_biology.toml' and add a
   constant beam matching the one used in the excellent paper of R. Labarbe.

### [bio_param]

Contains the main biological parameters (only pH for now)

### [concentrations.fixed]

Contains the concentrations of species supposed to be constant.
=> If a species is recorded in this section, its concentration will be forced constant.
For example:

```toml
[concentrations.fixed] # Unit is [mol]/[l]
H2O = 55
catalase = 8.0e-8
```

### [concentrations.initial]

Contains the initial concentrations of the species.
=> If a species is not recorded in this section, its initial concentration will be
assumed as 0.0 mol/l.
For example:

```toml
[concentrations.initial] # Unit is [mol]/[l]
O2 = 50e-6
```

### [reactions.radiolytic]

Contains the G-values of the species created by the radiation
source. The unit is [radical/100eV/incident particle].
For example:

```toml
[reactions.radiolytic]
e_aq = 2.8
OH_r = 2.8
```

### [[reactions.acid_base]]

Contains all the acid/base reactions involved in the model. Each
entry must defined 3 fields: "acid", "base" & "pKa".
For example:

```toml
[[reactions.acid_base]]
acid = "OH_r"
base = "O_r_minus"
pKa = 11.9

```

ps: note the double "[[" and "]]", mandatory for each entry!

### [[reactions.k_reaction]]

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
Reactants and products must only be separated by '->' characters.

Blanks/spaces for separation are allowed but not mandatory.

ps: note the double "[[" and "]]", mandatory for each entry!

### [[reactions.enzymatic]]

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
must be specified in the operator arrow (-- 'enzyme name' >>).

ps: note the double "[[" and "]]", mandatory for each entry!

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

ps: note the double "[[" and "]]", mandatory for each entry!

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

G-values are extrapolated according to the work of [D.
Boscolo](https://doi.org/10.3390/ijms21020424)  
Therefore, LET values must be in the range: [0.13, 150] keV/µm

## Ways of Running Simulation

The module currently supports 2 ways of running simulations:

- run in one single block

   ```python
   import radiopyo as rp 

   # Simulation handler
   sim = rp.UnitCell.from_toml_file(rp.LABARBE_CONFIG_FILE)
   res = sim.run(t_span=[1e-9, 10]) #Run simulation from 1e-9 s to 10 s.

   ```

- run in several chunks (useful for pulsed beams)

   ```python
   import radiopyo as rp 

   # Beam Definition
   peak_dose_rate  = 1000 # Gy/s
   max_dose        = 10   # Gy
   period          = 1e-3 # seconds
   on_time         = 5e-6 # seconds
   beam = rp.PulsedBeam.from_peak_dose_rate(peak_dose_rate, period, on_time, max_dose)
   
   # Simulation handler
   sim = rp.UnitCell.from_toml(rp.LABARBE_CONFIG_FILE)
   sim.set_beam(beam)
   
   # Better to adapt simulation settings during ON and OFF time.
   res = sim.prepare_chunked_run(t_span=[1e-9, 10],
                                 max_step_size_on=1e-8,  #seconds
                                 max_step_size_off=0.01, #seconds
                                 ).run()

   ```

## Acknowledgement

- This module was originally inspired by R. Labarbe
  [validated](https://doi.org/10.1016/j.radonc.2020.06.001) MatLab code.
