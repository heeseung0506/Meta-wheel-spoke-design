# Meta-Wheel Spoke Design (MATLAB)

This repository contains MATLAB scripts developed for the **3D spoke design of a meta-wheel** with 72 spoke pairs.  
The code generates the geometry, mesh, and finite element input files for penetration-preventing contact analysis.

## Features
- **Spoke geometry generation**: Defines backbone curve, thickness distribution, and 3D extrusion
- **Mesh generation**: Creates structured hexahedral mesh (nodes and C3D8R elements)
- **Export**: Writes ABAQUS `.inp` file for finite element simulations
- **Visualization**: Plots 3D spoke mesh for inspection and verification

## Requirements
- MATLAB R2023b (tested)
- No additional toolboxes required (base MATLAB is sufficient)

## Quick Start
1. Clone or download this repository:
   ```bash
   git clone https://github.com/heeseung0506/Meta-wheel-spoke-design.git

