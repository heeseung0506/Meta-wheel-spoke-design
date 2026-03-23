# Meta-Wheel Spoke Design (MATLAB)

> MATLAB scripts for parametric design, FE mesh generation, stiffness analysis, and Kriging surrogate modeling of a **meta-wheel** with curved spoke geometry.

---

## 📁 Repository Structure

| File | Description |
|------|-------------|
| `discrete_curved_spoke_design.m` | Generates 3D spoke geometry and exports ABAQUS `.inp` file for 72-spoke circular array |
| `Kriging_analytical.m` | Builds Kriging (GPR) surrogate models and plots 2-parameter design maps (K_Y, K_X, stiffness ratio) |
| `alpha_multi_axial_stiffness.m` | Parses FE results by helix angle α (r = 1, 3, 4, 6) and plots axial/lateral stiffness curves |
| `beta_multi_axial_stiffness.m` | Parses FE results by bias angle β (0°, 10°, 20°, 30°) and plots force–displacement curves |
| `radial_force_distribution.m` | Plots polar radial force distribution for continuous and discrete spoke configurations |
| `Inclined_Cantilever_Beam_Deflection_Visualization.m` | Analytical visualization of inclined cantilever beam deflection under vertical load |

---

## 🔩 Design Parameters

| Symbol | Description | Range |
|--------|-------------|-------|
| α (alpha) | Helix angle of spoke | 4.8° – 26.7° |
| β (beta) | Bias angle of spoke | 0° – 30° |
| r | Apex height of backbone curve | 1 – 6 mm |
| L | Vertical height of spoke | 75 mm |
| D | Depth of spoke | 60 mm |

---

## ✨ Features

- **Spoke geometry generation** — Defines backbone curve, thickness distribution, and 3D hexahedral mesh
- **ABAQUS export** — Writes `.inp` file for 72-spoke circular array (C3D8R elements)
- **Multi-axial stiffness analysis** — Extracts K_Y (axial) and K_X (lateral) stiffness from FE simulation data
- **Kriging surrogate modeling** — Fits GPR models (Matérn 5/2 kernel) over (α, β) parameter space
- **Design maps** — Contour plots of K_Y, K_X, and log₁₀(K_Y/K_X) for rapid design exploration
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
% Output: Figure 7 (a) K_Y map  (b) K_X map  (c) Stiffness ratio map
```

---

## 📊 Output Examples

### Kriging Design Maps (Figure 7)
| Panel | Variable | Description |
|-------|----------|-------------|
| (a) | K_Y | Axial stiffness [N/mm] — cool colormap |
| (b) | K_X | Lateral stiffness [N/mm] — hot colormap |
| (c) | log₁₀(K_Y / K_X) | Stiffness anisotropy ratio — RdYlGn colormap |

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
- **Hyperparameter optimization**: Bayesian optimization (40 evaluations)
- **Mesh type**: Structured hexahedral (C3D8R) for ABAQUS/Explicit
- **Stiffness extraction**: Initial slope (linear regression on first 30% of F-δ curve)

---

## 📬 Contact

**Heeseung** — [@heeseung0506](https://github.com/heeseung0506)

Feel free to open an issue if you have questions or suggestions.
