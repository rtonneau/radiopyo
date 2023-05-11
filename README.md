# radiopyo: Python Module dedicated to Radiobiology

![version](https://img.shields.io/badge/version-0.1.0-blue)
![License](https://img.shields.io/badge/License-LGPL3-radiopyo.svg)
![maintained](https://img.shields.io/badge/Maintained%3F-yes-green.svg)
![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v1.json)

This python3 module allows you to perform radiobiology related simulations. The embedded
science (all Physics, Chemistry & Biology) has the purpose to solve systems of Ordinary
Differential Equations (a.k.a. ODE) based on equations provided via an input file.

## Getting Started

This documentation assumes you already have or know how to install a new python3
environment on your local machine. If it is not yest the suggest, I suggest to check
[this](https://www.anaconda.com/) out.

Currently the main dependencies of this modules are:

- scipy
- numpy
- pandas
- pandas-stubs

## Installation

1. Run:

   ```bash
   pip install git+https://gitlab.unamur.be/rtonneau/radiopyo.git
   ```

2. You are good to go.

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
                                 max_step_size_on=1e-8,  #seconds
                                 max_step_size_off=0.01, #seconds
                                 ).run()

   df = res.to_pandas() # pandas.DataFrame
   
   # Time integrate species concentrations.
   s = res.integrate_species() # pandas.Series
   ```

### Configuring Simulation

For the moment, the only ways to load simulation configurations is via
[ron](https://github.com/ron-rs/ron) or [toml](https://toml.io/en/) files.
TOML files are recommended.
The TOML configuration file consists of several sections:

#### -> bio_param

Contains the main biological parameters (pH)

#### -> concentrations.fixed

Contains the concentrations of species supposed to be constant.
=> If a species is recorded in this section, its concentration will be forced constant.
For example:
```toml
[concentrations.fixed] # Unit is [mol]/[l]
H2O = 55
catalase = 8.0e-8
```

#### -> concentrations.initial

Contains the initial concentrations of the species.
=> If a species is not recorded in this section, its initial concentration will be
assumed as 0.0 mol/l.
For example:
```toml
[concentrations.initial] # Unit is [mol]/[l]
O2 = 50e-6
```

#### -> reactions.radiolytic

Contains the G values of the species created by the radiation
source. The unit is [radical/100eV/incident particle].
For example:
```toml
[reactions.radiolytic]
e_aq = 2.8
OH_r = 2.8
```

#### -> reactions.acid_base

Contains all the acid/base reactions involved in the model. Each
entry must defined 3 fields: "acid", "base" & "pKa".
For example:

```toml
[[reactions.acid_base]]
acid: "OH_r"
base: "O_r_minus"
pKa: 11.9

```

#### -> reactions.k_reaction

Contains all the simple chemical reactions involved in the model.
Each list entry must defined 3 fields: "reactants", "products" & "k_value".
For example:

```toml
[[reactions.k_reaction]]  # 25) O_r_minus + O2_r_minus -> 2 OH_minus + O2
reactants = ["O_r_minus", "O2_r_minus",]
products = ["OH_minus", "OH_minus"," O2",]
k_value = 6.0e8
```

#### -> enzymatic_reactions

Contains all the Enzymatic reactions involved in the model. The
enzymatic reactions are modelled with the Michaelis-Menten kinetic. Each list
entry must defined 5 fields: "enzyme", "substrate", "products", "k_value" & "k_micha".
For example:

```toml
[[reactions.enzymatic]]
enzyme =" catalase"
substrate = ["H2O2", " H2O2",]
products = ["O2", "H2O", "H2O",]
k_value = 6.62e7 #  [1/s]
k_micha = 1.1 #  [mol/l]

```

### Ways of Running Simulation

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

### Acknowledgement

- This module was originally inspired by R. Labarbe
  [validated](https://doi.org/10.1016/j.radonc.2020.06.001) MatLab code.