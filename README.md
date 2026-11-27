
# Investigation of True and Spurious Eigenvalues and Animation of Membrane Eigenmodes

## Authors

- Jeng-Tzong Chen
- Yi-Hsuan Yu
- Shing-Kai Kao
- Ying-Te Lee

Department of Harbor and River Engineering, National Taiwan Ocean University, Taiwan

## Overview

This repository provides the source code and numerical data associated with the manuscript entitled **"Investigation of True and Spurious Eigenvalues and Animation of Membrane Eigenmodes"**, submitted to the *Journal of Vibration Engineering & Technologies* (JVET).

The present study employs the dual boundary element method (BEM) combined with singular value decomposition (SVD) to investigate membrane eigenvalue problems governed by the two-dimensional Helmholtz equation under Dirichlet boundary conditions.

## Research Scope

- Eigenvalue search
- Identification of true and spurious eigenvalues
- Eigenvalue multiplicity determination
- Boundary discretization and convergence analysis
- Vibration-mode reconstruction
- Modal independence and orthogonality verification
- Taking-free-body analysis
- Modal animation

## Repository Structure

| Directory / File | Description | Manuscript Section |
|---|---|---|
| `Circle DIRICHLET` | Circular membrane | Section 3.1 |
| `Annulus DIRICHLET` | Annular membrane | Section 3.2 |
| `(3,1)-gon DIRICHLET` | Regular triangular membrane | Section 3.3 |
| `(5,2)-gon DIRICHLET` | Five-point star-shaped membrane | Section 3.4 |
| `(7,3)-gon DIRICHLET` | Seven-point star-shaped membrane | Section 3.5 |
| `L_shape DIRICHLET` | L-shaped membrane | Section 3.6 |
| `Taking free body` | Taking-free-body analysis | Section 3.7 |
| `(5,2)-gon (Circular Approximation) DIRICHLET` | Circular approximation of the nodal line | Section 3.7 |
| `(5,2)-gon (Nodal Line) DIRICHLET` | Numerically extracted nodal-line partitioning | Section 3.7 |
| `FEM_RegularPolygon_Eigen.m` | FEM-based numerical validation | Sections 3.4 and 3.5 |

## Programming Languages

- Fortran
- MATLAB

## File Types

- `.FOR`: Fortran source code
- `.m`: MATLAB scripts
- `.DAT`: Numerical input and output data
- `.TXT`: Supplementary numerical results

## General Workflow

1. Select the directory corresponding to the desired numerical example.
2. Prepare the boundary geometry and discretization using the corresponding MATLAB scripts.
3. Compile and execute the corresponding Fortran source code for eigenvalue analysis.
4. Examine the singular-value curves using the MATLAB plotting scripts.
5. Select the eigenvalue of interest and perform the corresponding vibration-mode calculation.
6. Visualize the boundary modes and interior vibration modes using MATLAB.

The specific input parameters and execution procedures may vary among numerical examples.

The eigenvalue search generally employs a uniform wavenumber increment of 0.001, as described in the manuscript.

## Notes

- Preserve the original directory structure and file names when running the programs.
- Numerical data files are included to facilitate examination of the computational results.
- The accompanying modal animations are described in the manuscript as supplementary online resources.

## Corresponding Author

Prof. Jeng-Tzong Chen  
National Taiwan Ocean University  
Email: jtchen@mail.ntou.edu.tw
