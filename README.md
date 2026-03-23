# Meta-Wheel Spoke Design (MATLAB)

> MATLAB scripts for parametric design, FE mesh generation, stiffness analysis, and Kriging surrogate modeling of a **meta-wheel** with curved spoke geometry.

---

## 📁 Repository Structure

| File | Description |
|------|-------------|
| `discrete_curved_spoke_design.m` | Generates 3D spoke geometry using a cosine-based backbone profile, creates structured hexahedral mesh (C3D8R) for 72 spokes in a circular array, and exports an ABAQUS-compatible `.inp` file |
| `Kriging_analytical.m` | Fits Gaussian Process Regression (Kriging) surrogate models to FE stiffness data over (α, β) parameter space and generates two-parameter design maps for K_X, K_Y, and their ratio |
| `alpha_multi_axial_stiffness.m` | Parses FE force-displacement data from Excel and computes axial/lateral stiffness for varying helix angle α (r = 1, 3, 4, 6 mm), extracting initial slope stiffness from the first 30% of the F-δ curve |
| `beta_multi_axial_stiffness.m` | Parses FE force-displacement data and evaluates multi-axial stiffness for varying bias angle β (0°–30°), plotting vertical, longitudinal, and lateral F-δ curves with optional displacement markers at U = 5 mm |
| `radial_force_distribution.m` | Loads radial force data from CSV files for continuous and discrete spoke configurations, and generates polar plots comparing force distribution across spoke types and β angles |
| `Inclined_Cantilever_Beam_Deflection_Visualization.m` | Analytical visualization of inclined cantilever beam deflection under vertical load |
| `Cantilever_Longitudinal,Lateral.m` | Implements inclined cantilever beam models for both longitudinal (Eq. 26: δ = F·cosβ·L³/3EI) and lateral (Eq. 30: δ = F·sinβ·L³/3EI) directions in a single 1×2 subplot, showing stiffness increase and compliance increase with β respectively |
| `Longitudinal Stiff Eq.(10).m` | Computes normalized effective longitudinal stiffness k_eff(α) using the axial-bending projection model (Eq. 10), marking reference design points at α = 4.8°, 14.1°, 18.5°, 26.7° and quantifying the 17.7% stiffness drop across the studied curvature range |

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

## ✨ Features

- **Spoke geometry generation** — Defines cosine-based backbone curve, thickness distribution, and 3D hexahedral mesh
- **ABAQUS export** — Writes `.inp` file for 72-spoke circular array (C3D8R elements)
- **Multi-axial stiffness analysis** — Extracts K_Y (lateral) and K_X (longitudinal) stiffness from FE simulation data
- **Kriging surrogate modeling** — Fits GPR models (Matérn 5/2 kernel) over (α, β) parameter space
- **Design maps** — Contour plots of K_Y, K_X, and log₁₀(K_Y/K_X) for rapid inverse design
- **Analytical beam models** — Closed-form stiffness expressions (Eq. 10, 26, 30) for longitudinal and lateral directions
- **Radial force visualization** — Polar plots comparing continuous and discrete spoke configurations

---

## ⚙️ Requirements

| Script | MATLAB Version | Required Toolbox |
|--------|---------------|-----------------|
| `discrete_curved_spoke_design.m` | R2021b+ | None (base MATLAB) |
| `Kriging_analytical.m` | R2021b+ | **Statistics and Machine Learning Toolbox** (`fitrgp`) |
| `alpha_multi_axial_stiffness.m` | R2021b+ | None |
| `beta_multi_axial_stiffness.m` | R2021b+ | None |
| `radial_force_distribution.m` | R2021b+ | None |
| `Inclined_Cantilever_Beam_Deflection_Visualization.m` | R2021b+ | None |
| `Cantilever_Longitudinal,Lateral.m` | R2021b+ | None |
| `Longitudinal Stiff Eq.(10).m` | R2021b+ | None |

> Tested on **MATLAB R2023b**.

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

**3. Run stiffness analysis (alpha-sweep)**
```matlab
% Place stiffness.xlsx in the working directory
run('alpha_multi_axial_stiffness.m')
```

**4. Run Kriging design maps**
```matlab
% Requires Statistics and Machine Learning Toolbox
run('Kriging_analytical.m')
% Output: Figure 7 — (a) K_X map  (b) K_Y map  (c) Stiffness ratio map
```

**5. Run analytical beam models**
```matlab
run('Longitudinal Stiff Eq.(10).m')       % Figure 5b  — k_eff vs α (Eq. 10)
run('Cantilever_Longitudinal,Lateral.m')  % Figure 10b & 11b — longitudinal & lateral vs β
```

---

## 📊 Output Examples

### Kriging Design Maps (Figure 7)
| Panel | Variable | Description |
|-------|----------|-------------|
| (a) | K_X | Longitudinal stiffness [N/mm] |
| (b) | K_Y | Lateral stiffness [N/mm] |
| (c) | log₁₀(K_Y / K_X) | Stiffness anisotropy ratio — RdYlGn colormap |

### Analytical Model Summary

| Script | Equation | Key Result |
|--------|----------|------------|
| `Longitudinal Stiff Eq.(10).m` | Eq. (10): k_eff = k_ax·cos²α + k_b·sin²α | 17.7% stiffness drop from α = 4.8° to 26.7° |
| `Cantilever_Longitudinal,Lateral.m` | Eq. (26): δ = F·cosβ·L³/3EI | Longitudinal stiffness increases with β |
| `Cantilever_Longitudinal,Lateral.m` | Eq. (30): δ = F·sinβ·L³/3EI | Lateral compliance increases with β |

### Stiffness Data (FE simulation)

**K_Y [N/mm]** — varies with α and β:

| α \ β | 0° | 10° | 20° | 25° | 30° |
|--------|-----|-----|-----|-----|-----|
| 4.8° | 46.15 | 44.88 | 41.21 | 38.58 | 35.49 |
| 14.1° | 24.76 | 17.48 | 14.98 | 14.18 | 13.38 |
| 18.5° | 21.62 | 16.43 | 15.54 | 14.99 | 14.42 |
| 22.7° | 19.53 | 18.09 | 16.69 | 15.80 | 14.82 |
| 26.7° | 17.24 | 16.76 | 15.42 | 14.64 | 13.62 |

---

## 📂 Required Data Files

The following external data files are needed for stiffness analysis scripts:

| File | Used by |
|------|---------|
| `stiffness.xlsx` | `alpha_multi_axial_stiffness.m` |
| `beta_stiffness.xlsx` | `beta_multi_axial_stiffness.m` |
| `contin.csv`, `0degree.csv`, `10degree.csv`, `20degree.csv`, `2DTWEEL2.csv` | `radial_force_distribution.m` |

> ⚠️ Data files are not included in this repository (FE simulation outputs). Contact the author if needed.

---

## 📌 Method Summary

- **Surrogate model**: Gaussian Process Regression (Kriging) with Matérn 5/2 kernel
- **Hyperparameter optimization**: Bayesian optimization (40 evaluations), 5-fold CV for K_Y, 4-fold CV for K_X
- **Mesh type**: Structured hexahedral (C3D8R) for ABAQUS/Explicit
- **Stiffness extraction**: Initial slope (linear regression on first 30% of F-δ curve)
- **Analytical model**: Euler-Bernoulli cantilever beam with axial-bending projection (Eqs. 10, 26, 30)

---

## 📄 Related Publication

This repository accompanies the following manuscript:

> H.Han, "Directional Stiffness Decoupling in Meta-Wheels via Three-Dimensional Discrete Curved Spokes," *under review*, 2026.

> ⚠️ MATLAB codes will be made publicly available upon acceptance of the manuscript, or provided upon reasonable request to the corresponding author.

---

## 📬 Contact

**Heeseung** — [@heeseung0506](https://github.com/heeseung0506)

Feel free to open an issue if you have questions or suggestions.
