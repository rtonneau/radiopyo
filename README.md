# radiopyo: Python Module dedicated to Radiobiology

![version](https://img.shields.io/badge/version-0.1.0-blue)
[![npm license](https://img.shields.io/npm/l/badge-maker.svg)](https://npmjs.org/package/badge-maker)
![maintained](https://img.shields.io/badge/Maintained%3F-yes-green.svg)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v1.json)](https://github.com/charliermarsh/ruff)
[![image](https://img.shields.io/pypi/l/ruff.svg)](https://pypi.python.org/pypi/ruff)

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

   # Beam Definition
   dose_rate  = 1000 # Gy/s
   max_dose = 10 # Gy
   beam = rp.ConstantBeam(dose_rate, max_dose)

   # Simulation handler
   sim = rp.UnitCell.from_ron_file(rp.LABARBE_CONFIG_FILE, beam)
   res = sim.run(t_span=[1e-9, 10]) #Run simulation from 1e-9 s to 10 s.

   # Extract results to a nice pandas DataFrame.
   df = res.to_pandas() 
   ```

### Configuring Simulation

For the moment, the only way to load simulation configuration is via
[ron](https://github.com/ron-rs/ron) files. The configuration file consists of several
sections: 

#### -> bio_param

Dictionary like structure holding the main biological parameters (pH)

#### -> fixed_concentrations

Dictionary like structure holding the concentrations of species supposed to be constant.
=> If a species is recorded in this section, its concentration will be forced constant.

#### -> initial_concentrations

Dictionary like structure holding the initial concentrations of the species.
=> If a species is not recorded in this section, its initial concentration will be
assumed as 0.0 mol/l.

#### -> radiolytic

Dictionary like structure holding the G values of the species created by the radiation
source. The unit is [radical/100eV/incident particle].

#### -> acid_base

List like structure holding all the acid/base reactions involved in the model. Each list
entry must defined 3 fields: "acid", "base" & "pKa".
For example:

```yaml
(
   acid: "OH_r",
   base: "O_r_minus",
   pKa: 11.9,
),
```

#### -> k_reactions

List like structure holding all the simple chemical reactions involved in the model.
Each list entry must defined 3 fields: "reactants", "products" & "k_value".
For example:

```yaml
//O_r_minus + O2_r_minus -> 2 OH_minus + O2
(
   reactants: ["O_r_minus", "O2_r_minus"],
   products: ["OH_minus", "OH_minus", "O2"],
   k_value: 6.0e8
),
```

#### -> enzymatic_reactions

List like structure holding all the Enzymatic reactions involved in the model. The
enzymatic reactions are modelled with the Michaelis-Menten kinetic. Each list
entry must defined 5 fields: "enzyme", "substrate", "products", "k_value" & "k_micha".
For example:

```yaml
(
   enzyme: "catalase",
   substrate: ["H2O2", "H2O2"],
   products: ["O2", "H2O", "H2O"],
   k_value: 6.62e7, // [1/s]
   k_micha: 1.1, // [mol/l]
),
```

### Ways of Running Simulation

The module currently supports 2 ways of running simulations:

- run in one single block

   ```python
   import radiopyo as rp 

   # Beam Definition
   dose_rate  = 1000 # Gy/s
   max_dose = 10 # Gy
   beam = rp.ConstantBeam(dose_rate, max_dose)

   # Simulation handler
   sim = rp.UnitCell.from_ron_file(rp.LABARBE_CONFIG_FILE, beam)
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
   sim = rp.UnitCell.from_ron_file(rp.LABARBE_CONFIG_FILE, beam)

   chunks = sim.beam.make_time_intervals(start=1e-9,stop=500)
   #chunks[:, 1] contains the max_step variable of ODE solver
   chunks[np.where(chunks[:,1]>0)[0], 1] = sim.beam.timer.on_time / 1e4
   chunks[np.where(chunks[:,1]==0)[0], 1] = (sim.beam.timer.period-sim.beam.timer.on_time) / 1e2
   if chunks[-1,0] > sim.beam.time_at_max_dose:
      chunks[-2:,1] = 0.01 
   to_run = sim.prepare_chunked_run(chunks[:,0], max_step=chunks[:,1], use_jac=True,
   method="LSODA", atol=1e-10)
   res = ro_tun.run_all()

   ```

### Acknowledgement

- This module was originally inspired by R. Labarbe
  [validated](https://doi.org/10.1016/j.radonc.2020.06.001) MatLab code.