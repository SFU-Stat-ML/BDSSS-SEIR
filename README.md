# Beta-Dirichlet Switching State-Space SEIR Model (BDSSS-SEIR)
**Goal**: We propose a Beta-Dirichlet switching state-space transmission model to track underlying dynamics of disease and evaluate the effectiveness of interventions simultaneously.

## Features
This repository contains R code for simulation studies and real data analysis of BDSSS-SEIR model, described as follows:
- "Code" contains source R code for particle Gibbs samplers in two-regime, three-regime and four-regime settings.
  - "Two Regime":
    - Users can run "SimulationStudy_TwoRegimes.R" to simulate two-regime data from the BDSSS-SEIR model and implement "particleGibbs_TwoRegimes.R" for inference.
    - Users can run "BCData_TwoRegime.R" to read in BC COVID-19 observed infectious proportions and implement "particleGibbs_TwoRegimes.R" for inference.
  - "Three Regime":
    - Users can run "SimulationStudy_ThreeRegimes.R" to simulate three-regime data from the BDSSS-SEIR model and implement "particleGibbs_ThreeRegimes.R" for inference.
    - Users can run "BCData_ThreeRegime.R" to read in BC COVID-19 observed infectious proportions and implement "particleGibbs_ThreeRegimes.R" for inference.
  - "Four Regime":
    - Users can run "BCData_FourRegime.R" to read in BC COVID-19 observed infectious proportions and implement "particleGibbs_FourRegimes.R" for inference.
- "Data" contains simulated data in two-regime and three-regime settings, as well as the real data retreived from [here](https://docs.google.com/spreadsheets/d/1KvX2bNs4hUYGY8Kk47c4SmFVHrKT_0vEYY0NqKQShZs/edit#gid=0).

