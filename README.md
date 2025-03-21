# Bismuth-207 Purity Monitor Simulation & Analysis
This repository contains a set of tools for simulating, processing, and analyzing the decay spectrum of a Bi-207 radioactive source in a dual purity monitor (PM) setup. The setup models two concentric rings (inner and outer) and is used to understand signal behavior captured by an oscilloscope.

---

###  Program List and Descriptions

#### Simulation

- **`bi207PMtime.f`**
  - Original Fortran code by Francesco.
  - Simulates the decay of Bi-207 and records signal hits in two rings (inner and outer) of the purity monitors.
  - **Output**: `bi207stream.txt`

- **`bi207PMtime.py`**
  - Python translation of `bi207PMtime.f` (used because the Fortran script couldn't be executed on some systems).
  - Intended to replicate `bi207PMtime.f`, but may have slight differences (to be validated).
  - **Output**: `bi207stream.txt`

#### Scope Emulation (Oscilloscope-like Behavior)

- **`bi207PMscope.f`**
  - Original Fortran code by Francesco.
  - Reads `bi207stream.txt` and simulates the logic of oscilloscope data acquisition.
  - Produces signal readouts per ring and drift length.
  - **Outputs**:
    - `inner_anode_long.txt`
    - `inner_anode_short.txt`
    - `outer_anode_long.txt`
    - `outer_anode_short.txt`
    - Combined: `bi207spectra.txt`

- **`bi207PMscope.py`**
  - Python translation of `bi207PMscope.f`.
  - Emulates the oscilloscope behavior described above in a more accessible language.

#### Analysis & Visualization

- **`calibration.py`**
  - Plots calibration histograms and fits Gaussian curves to determine uncertainties.
  - Used for test pulse calibration studies.

- **`plot_spectra.py`**
  - Plots energy spectra from `bi207spectra.txt`.
  - Compares simulated spectra with experimental data.

- **`plot_trace.py`**
  - Optional utility to visualize raw signal traces.
  - Reads from oscilloscope-style outputs:
    - `inner_anode_long.txt`
    - `inner_anode_short.txt`
    - `outer_anode_long.txt`
    - `outer_anode_short.txt`

---

###  File Associations

| Program File        | Inputs                     | Outputs                         |
|---------------------|----------------------------|----------------------------------|
| `bi207PMtime.f/py`  | Physical parameters        | `bi207stream.txt`               |
| `bi207PMscope.f/py` | `bi207stream.txt`          | Ring-wise `.txt` files + `bi207spectra.txt` |
| `calibration.py`    | Calibration pulses         | Gaussian fits / uncertainties    |
| `plot_spectra.py`   | `bi207spectra.txt`         | Spectral comparison plots        |
| `plot_trace.py`     | Ring `.txt` trace files    | Trace plots                      |

---


