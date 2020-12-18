# TransPore Version 1.0  
Date: 01/25/2017

A FEM-FDM solver for two-phase, multicomponent transport in porous media

## GENERAL USAGE NOTES
=====================

Requirements: MATLAB 2017a or higher, Windows/Linux/MacOS

E-mail: sdutta.math@gmail.com, daripa@math.tamu.edu

Copyright 2010-2018 TransPore developers and contributors. All rights reserved.

References: 

Daripa, P. & Dutta, S. (2017) Modeling and simulation of surfactant–polymer flooding using a new hybrid method. J. Comput. Phys., 335, 249–282.
https://doi.org/10.1016/j.jcp.2017.01.038

Daripa, P. & Dutta, S. (2017) Convergence analysis of a characteristics-based hybrid method for multicomponent transport in porous media, arXiv:1707.00035v1 [math.NA], 1–30.

Funding for this research was provided by:
Qatar National Research Fund (08-777-1-141),
National Science Foundation (DMS-1522782)


## Primary SOURCE FILE
===================

### Master_surf_grid.m
Variables and Data Structure:
nsim - number of different flooding simulations

sizeofgrid - Nx x Ny grid sizes for each simulation

c0iter, g0iter - nsim arrays of concentrations of components 1 & 2 respectively in the injected fluid for each simulation

f - source term for the elliptic problem. Nx x Ny matrix with non-zero intensities at injection and production wells

KK - Nx x Ny matrix with absolute permeability values for the domain

UU, CC, GG - Nx x Ny matrices for space-time values of wetting phase saturation, components 1 & 2 concentrations respectively 

miuw, miuo - wetting and non-wetting fluid base viscosities respectively

swr0, sor0 - wetting and non-wetting phase initial residual saturations respectively

sigma - Nx x Ny matrix for interfacial tension values over the domain

miua - Nx x Ny matrix for aqueous phase saturation values over the domain

lamba_a, lambda_o - Nx x Ny matrices for wetting phase and non-wetting phase mobility values over the domain 

u,v - Nx x Ny matrices for x-direction and y-direction total velocity values over the domain. Note: These are obtained by solving the global pressure equation and are different from phase velocities


## Secondary SOURCE FILES
=====================

KKdef() - function implementing different types of homogeneous, heterogeneous, stochastic, piecewise constant absolute permeability profiles

s0c0() - function implementing initial configurations and injection profiles for each simulation.

compvis(), compres(), compmob() - functions to update phase viscosities, phase residual saturations and phase mobility values with the evolution of state variables

setGrid(), setRightHand(), setA(), setB(), getu(), get_vn() - functions implementing various parts of the elliptic solver for the global pressure and total velocities

nmmoc_surf_mod_neumann() - function implementing the MMOC-FD procedure for solving the component transport equations 
  
  



