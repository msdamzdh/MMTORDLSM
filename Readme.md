# Multi-Material Topology Optimization Using Reaction-Diffusion Level Set Method

This project implements a multi-material topology optimization algorithm using a reaction-diffusion level set method for structural design problems.

## Overview

The code performs topology optimization on a cantilever beam example, optimizing the distribution of multiple materials to minimize compliance while satisfying volume constraints. It uses a level set method combined with a reaction-diffusion equation to evolve the material boundaries.

## Key Features

- Multi-material topology optimization 
- Reaction-diffusion level set method
- Finite element analysis
- Augmented Lagrangian optimization
- Visualization of optimization iterations

## Files

- `CanteliverExample.m`: Main script to set up and run the cantilever beam optimization example
- `TO_MMRDLSM.m`: Core topology optimization implementation

## Usage

1. Ensure MATLAB is installed
2. Open and run `CanteliverExample.m`
3. Adjust parameters in the script to modify the optimization problem

## Parameters

Key parameters that can be adjusted include:

- Material properties (`MP` struct)
- Geometry and boundary conditions (`GP` struct) 
- Optimization parameters (`OP` struct)

## Output

The code generates:

- Plot of the optimized topology at each iteration
- Video file showing the optimization progression
- Final compliance and volume fraction values

## Dependencies

- MATLAB (developed with R2021a, but should work with other recent versions)

## Result


https://github.com/user-attachments/assets/d697e76f-bc60-46d7-ad16-1617c87757b2

