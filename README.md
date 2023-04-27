# radiopyo: Python Module dedicated to Radiobiology

This python3 module allows you to perform radiobiology related simulations. The embedded
science (all Physics, Chemistry & Biology) has the purpose to solve systems of Ordinary
Differential Equations (a.k.a. ODE) based on equations provided via an input file.


## Getting Started

## Installation

1. Clone or download the repository to your local machine.

    ```bash
    git clone https://github.com/SpringNO26/radiopyo.git
    ```

2. Enter the working directory.

   ```bash
   cd radiopyo
   ```

3. Install the required libraries.

   ```bash
   pip install -r requirements.txt
   # Or if you use mamba 
   mamba install -r requirements.txt
   ```

4. You are good to go.

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

### Ways of Running Simulation

## Project Organization

### Acknowledgement

* This module was originally inspired by R. Labarbe
  [validated](https://doi.org/10.1016/j.radonc.2020.06.001) MatLab code.