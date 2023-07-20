# Beta-Dirichlet Switching State-Space SEIR Model (BDSSS-SEIR)
### Authors: Jingxue Feng (jingxuef@sfu.ca), Liangliang Wang (lwa68@sfu.ca)
## Goal 
We propose a Beta-Dirichlet switching state-space transmission model to track underlying dynamics of disease and evaluate the effectiveness of interventions simultaneously. Bayesian inference are performed using a particle MCMC algorithm.

## Model Description
```math
   \boldsymbol{\theta}_t | \boldsymbol{\theta}_{t-1},x_t \sim \text{Dirichlet}(\kappa  r(\boldsymbol{\theta}_{t-1}; \alpha, \beta, \gamma, f_{x_t})) \\
```
```math
   y_t |\boldsymbol{\theta}_t \sim \text{Beta}(\lambda p I_t, \lambda (1- p I_t))
```
- $y_t$ is the observed infectious proportion at time $t$ ;
- $x_t$ is a discrete latent variable describing the regime of model at time $t$;
- $p$ is the identification rate representing the fraction of infectious population identified by diagnosis;
- $\boldsymbol{\theta}_t = [S_t, E_t, I_t, R_t]^\top$ is the vector of the underlying proportion of the susceptible, exposed, infectious, and recovered population, where $S_t+E_t+I_t+R_t=1 \ \forall t$;
- $\alpha$, $\beta$ and $\gamma$ are the latency rate, transmission rate and recovery rate in the SEIR system respectively;
- $\lambda$ and $\kappa$ are the precision parameters controlling the variances (or randomness) for the observation and transition process;
- $f_{x_t}$ is the transmission rate modifier that changes with respect to the state variable $x_t$. 
- $\psi$ represents the set of model parameters to be estimated.
- $r(.)$ is the solution of the modified SEIR system starting at $\boldsymbol{\theta}_{t-1}$.
  
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
