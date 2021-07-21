# FV3 Examples

This directory contains Python (Jupyter) notebooks demonstrating basic FV3 capabilities, including characteristics of the solver, physics-dynamics coupling, output using FMS diag_manager, and basic analyses of the model output. The notebooks should all be viewable in any browser; you can also download any of them and use them in an up-to-date Python/Jupyter environment.

## 1D Notebooks
tp_core
: 1D advection operators in FV3. This is designed to be an *interactive* notebook for downloading and playing with the options, initial conditions, zooms, and so on.

## 2D Shallow-water Notebooks
RHwave
: Rossby-Haurwitz wave, a test of height-vorticity consistency

BTwave
: Barotropic instability, demonstrating wave breaking (cf. Galewsky et al. 2004; Scott, Harris, and Polvani 2016)

BLvortex
: Bates-Li forced polar vortex

## 3D Notebooks

BCwave
: Hydrostatic baroclinic wave, with and without moisture

TC
: Reed-Jablonowski TC tests, demonstrating the effect of advection schemes and numerical diffusion

DPsupercell
: Supercell on a doubly-periodic using Weisman (WK82) sounding and a straight-line hodograph. Demonstrates pressure partitioning between moist, dry, and nonhydrostatic contributions.
