<div align="center">

# 🏎 Meta-Wheel Spoke Design 🏎

**MATLAB scripts for parametric design, FE mesh generation, stiffness analysis,**  
**and Kriging surrogate modeling of a meta-wheel with 3D discrete curved spokes.**

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021b%2B-orange?logo=mathworks&logoColor=white)](https://www.mathworks.com/)
[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc/4.0/)
[![Status](https://img.shields.io/badge/Status-Under%20Review-yellow)]()
[![ABAQUS](https://img.shields.io/badge/FE%20Solver-ABAQUS-blue)]()

</div>

---

## 📌 Table of Contents
- [Overview](#-overview)
- [Repository Structure](#-repository-structure)
- [Design Parameters](#-design-parameters)
- [Requirements](#-requirements)
- [Quick Start](#-quick-start)
- [Output Examples](#-output-examples)
- [Required Data Files](#-required-data-files)
- [Method Summary](#-method-summary)
- [Related Publication](#-related-publication)
- [Contact](#-contact)

---

## 🔍 Overview

This repository provides MATLAB scripts developed for the study of **directional stiffness decoupling** in non-pneumatic meta-wheels. The design is governed by two independent geometric parameters:

| Parameter | Role |
|-----------|------|
| **α** — In-plane curvature angle | Controls geometric softening across all directions |
| **β** — Out-of-plane inclination angle | Enables decoupled tuning of longitudinal vs. lateral stiffness |

> **Key result:** Increasing β raises longitudinal stiffness by ~34% and reduces lateral stiffness by ~35%, while keeping vertical stiffness variation within 1.4%.

---

## 📁 Repository Structure

| File | Figure | Description |
|------|--------|-------------|
| `discrete_curved_spoke_design.m` | Fig. 2 | Generates 3D spoke geometry using a cosine-based backbone profile, creates structured hexahedral mesh (C3D8R) for 72 spokes in a circular array, and exports an ABAQUS-compatible `.inp` file |
| `alpha_multi_axial_stiffness.m` | Fig. 3, 4 | Parses FE force-displacement data from Excel and computes axial/lateral stiffness for varying helix angle α (r = 1, 3, 4, 6 mm), extracting initial slope stiffness from the first 30% of the F-δ curve |
| `beta_multi_axial_stiffness.m` | Fig. 6 | Parses FE force-displacement data and evaluates multi-axial stiffness for varying bias angle β (0°–30°), plotting vertical, longitudinal, and lateral F-δ curves with optional displacement markers at U = 5 mm |
| `Longitudinal Stiff Eq.(10).m` | Fig. 5b | Computes normalized effective longitudinal stiffness k_eff(α) using the axial-bending projection model (Eq. 10), marking reference design points at α = 4.8°, 14.1°, 18.5°, 26.7° and quantifying the 17.7% stiffness drop |
| `Kriging_analytical.m` | Fig. 7 | Fits Gaussian Process Regression (Kriging) surrogate models to FE stiffness data over (α, β) parameter space and generates two-parameter design maps for K_X, K_Y, and their ratio |
| `Cantilever_Longitudinal,Lateral.m` | Fig. 10b, 11b | Implements inclined cantilever beam models for both longitudinal (Eq. 26) and lateral (Eq. 30) directions in a single 1×2 subplot |
| `radial_force_distribution.m` | Fig. 9c | Loads radial force data from CSV files for continuous and discrete spoke configurations, and generates polar plots comparing force distribution across spoke types and β angles |
| `Inclined_Cantilever_Beam_Deflection_Visualization.m` | — | Analytical visualization of inclined cantilever beam deflection under vertical load |

---

## 🔩 Design Parameters

| Symbol | Description | Range |
|--------|-------------|-------|
| α (alpha) | In-plane curvature angle of spoke | 4.8° – 26.7° |
| β (beta) | Out-of-plane inclination angle of spoke | 0° – 30° |
| r | Apex height of backbone curve | 1 – 6 mm |
| L | Vertical height of spoke | 75 mm |
| D | Depth of spoke | 60 mm |

---

## ⚙️ Requirements

| Script | MATLAB | Toolbox |
|--------|--------|---------|
| `discrete_curved_spoke_design.m` | R2021b+ | None |
| `Kriging_analytical.m` | R2021b+ | ⚠️ **Statistics and Machine Learning Toolbox** |
| `alpha_multi_axial_stiffness.m` | R2021b+ | None |
| `beta_multi_axial_stiffness.m` | R2021b+ | None |
| `radial_force_distribution.m` | R2021b+ | None |
| `Longitudinal Stiff Eq.(10).m` | R2021b+ | None |
| `Cantilever_Longitudinal,Lateral.m` | R2021b+ | None |
| `Inclined_Cantilever_Beam_Deflection_Visualization.m` | R2021b+ | None |

> ✅ Tested on **MATLAB R2023b**

---

## 🚀 Quick Start

**1. Clone the repository**
```bash
git clone https://github.com/heeseung0506/Meta-wheel-spoke-design.git
cd Meta-wheel-spoke-design
```

**2. Generate ABAQUS input file**
```matlab
run('discrete_curved_spoke_design.m')
% Output: 72spokes_circular_array.inp
```

**3. Run stiffness analysis**
```matlab
run('alpha_multi_axial_stiffness.m')   % α-sweep (place stiffness.xlsx in folder)
run('beta_multi_axial_stiffness.m')    % β-sweep (place beta_stiffness.xlsx in folder)
```

**4. Run analytical beam models**
```matlab
run('Longitudinal Stiff Eq.(10).m')       % Figure 5b  — k_eff vs α
run('Cantilever_Longitudinal,Lateral.m')  % Figure 10b & 11b — stiffness vs β
```

**5. Run Kriging design maps**
```matlab
% ⚠️ Requires Statistics and Machine Learning Toolbox
run('Kriging_analytical.m')
% Output: Figure 7 — K_X map, K_Y map, Stiffness ratio map
```

---

## 📊 Output Examples

<details>
<summary><b>Kriging Design Maps (Figure 7)</b></summary>

| Panel | Variable | Description |
|-------|----------|-------------|
| (a) | K_X | Longitudinal stiffness [N/mm] |
| (b) | K_Y | Lateral stiffness [N/mm] |
| (c) | log₁₀(K_Y / K_X) | Stiffness anisotropy ratio — RdYlGn colormap |

</details>

<details>
<summary><b>Analytical Model Summary</b></summary>

| Script | Equation | Key Result |
|--------|----------|------------|
| `Longitudinal Stiff Eq.(10).m` | k_eff = k_ax·cos²α + k_b·sin²α | 17.7% stiffness drop: α = 4.8° → 26.7° |
| `Cantilever_Longitudinal,Lateral.m` | δ = F·cosβ·L³/3EI (Eq. 26) | Longitudinal stiffness increases with β |
| `Cantilever_Longitudinal,Lateral.m` | δ = F·sinβ·L³/3EI (Eq. 30) | Lateral compliance increases with β |

</details>

<details>
<summary><b>FE Stiffness Data — K_Y [N/mm]</b></summary>

| α \ β | 0° | 10° | 20° | 25° | 30° |
|--------|-----|-----|-----|-----|-----|
| 4.8° | 46.15 | 44.88 | 41.21 | 38.58 | 35.49 |
| 14.1° | 24.76 | 17.48 | 14.98 | 14.18 | 13.38 |
| 18.5° | 21.62 | 16.43 | 15.54 | 14.99 | 14.42 |
| 22.7° | 19.53 | 18.09 | 16.69 | 15.80 | 14.82 |
| 26.7° | 17.24 | 16.76 | 15.42 | 14.64 | 13.62 |

</details>

---

## 📂 Required Data Files

> ⚠️ Data files are **not included** in this repository (FE simulation outputs).  
> Contact the author if needed.

| File | Used by |
|------|---------|
| `stiffness.xlsx` | `alpha_multi_axial_stiffness.m` |
| `beta_stiffness.xlsx` | `beta_multi_axial_stiffness.m` |
| `contin.csv`, `0degree.csv`, `10degree.csv`, `20degree.csv`, `2DTWEEL2.csv` | `radial_force_distribution.m` |

---

## 📌 Method Summary

- 🔷 **Surrogate model** — Gaussian Process Regression (Kriging) with Matérn 5/2 kernel
- 🔷 **Hyperparameter optimization** — Bayesian optimization (40 evaluations), 5-fold CV for K_Y, 4-fold CV for K_X
- 🔷 **Mesh type** — Structured hexahedral (C3D8R) for ABAQUS/Standard
- 🔷 **Stiffness extraction** — Initial slope via linear regression on first 30% of F-δ curve
- 🔷 **Analytical model** — Euler-Bernoulli cantilever beam with axial-bending projection (Eqs. 10, 26, 30)

---

## 📄 Related Publication

> H.Han, "Directional Stiffness Decoupling in Meta-Wheels via Three-Dimensional Discrete Curved Spokes," *under review*, 2026.

> ⚠️ MATLAB codes will be made publicly available upon acceptance, or provided upon reasonable request to the corresponding author.

---

## 📬 Contact

<div align="center">

**Heeseung Han** — [@heeseung0506](https://github.com/heeseung0506)

*Feel free to open an issue if you have questions or suggestions.*

</div>
