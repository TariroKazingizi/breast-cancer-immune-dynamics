# MATLAB Simulations

This folder contains the MATLAB scripts used for the computational analysis of the breast cancer–immune dynamics research project.

The simulations investigate tumour growth, immune-cell responses, chemotherapy treatment strategies, stability, bifurcation behaviour, and tumour persistence.

## Simulations

### 1. Early vs Late Detection
**File:** `early_vs_late_detection.m`

Compares tumour dynamics under early and late detection scenarios and illustrates how the timing of intervention affects tumour progression.

### 2. Chemotherapy Comparison
**File:** `chemotherapy_comparison.m`

Compares constant and adaptive chemotherapy strategies. The simulation also examines immune-cell dynamics, drug concentration and the effect of varying the tumour growth rate.

### 3. Continuous vs Pulsed Chemotherapy
**File:** `continuous_vs_pulsed.m`

Compares continuous and pulsed chemotherapy treatment strategies using tumour dynamics, immune-cell dynamics, drug concentration and phase-plane analysis.

### 4. Stability and Coexistence
**File:** `stability_coexistence.m`

Investigates the stability of the tumour–immune system and the coexistence behaviour of tumour and immune-cell populations.

### 5. Hopf Bifurcation and Limit Cycle
**File:** `hopf_bifurcation_limit_cycle.m`

Explores Hopf bifurcation behaviour and the emergence of limit-cycle dynamics in the mathematical model.

### 6. Stability and Bifurcation Analysis
**File:** `bifurcation_eigenvalue_analysis.m`

Analyses equilibrium behaviour, Jacobian eigenvalues, stability and bifurcation patterns of the model.

### 7. Tumour Persistence vs Tumour-Free Dynamics
**File:** `tumor_persistence_vs_tumor_free.m`

Examines conditions associated with tumour persistence and tumour-free outcomes.

### 8. Coexistence and Phase-Plane Analysis
**File:** `breast_cancer_coexistence_short.m`

Investigates coexistence between tumour and immune populations using numerical simulation and phase-plane analysis.

### 9. Tumour Threshold Dynamics
**File:** `basin_of_attraction_tumor_threshold.m`

Investigates tumour threshold behaviour by examining how different initial tumour sizes influence the final tumour population.

## Computational Methods

The simulations use MATLAB numerical methods, including:

- Ordinary differential equation (ODE) solving
- `ode45` numerical integration
- Parameter variation
- Equilibrium analysis
- Jacobian and eigenvalue analysis
- Stability and bifurcation analysis
- Phase-plane analysis
- Chemotherapy treatment simulations

## Model Components

The mathematical model considers interactions between:

- Tumour cells
- Cytotoxic T lymphocytes (CTLs)
- Natural killer (NK) cells
- Gamma-delta T cells
- Chemotherapy drug concentration

## Reproducibility

The MATLAB scripts are provided to document and reproduce the computational simulations presented in the undergraduate research project.

Each script can be opened and executed in MATLAB independently.
