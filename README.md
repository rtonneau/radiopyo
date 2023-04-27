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
   from radiopyo import UnitCell, LABARBE_CONFIG_FILE

   # Beam Definition
   dose_rate  = 1000 # Gy/s
   max_dose = 10 # Gy
   beam = Constantbeam(dose_rate, max_dose)

   # Simulation handler
   sim = UnitCell.from_ron_file(LABARBE_CONFIG_FILE, beam)
   res = sim.run(t_span=[1e-9, 10]) #Run simulation from 1e-9 s to 10 s.

   # Extract results to a nice pandas DataFrame.
   df = res.to_pandas() 
   ```

### Configuring Simulation

### Running Simulation

## Project Organization

### Acknowledgement

* This module is inspired by R. Labarbe MatLab code and validated [here](https://doi.org/10.1016/j.radonc.2020.06.001)