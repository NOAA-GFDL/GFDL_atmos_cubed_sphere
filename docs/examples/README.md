# FV3 Examples

This directory contains Python (Jupyter) notebooks demonstrating basic FV3 capabilities, including characteristics of the solver, physics-dynamics coupling, output using FMS `diag_manager`, and basic analyses of the model output. The notebooks should all be viewable in any browser; you can also download any of them and use them in an up-to-date Python/Jupyter environment.

## 1D Cases
tp_core
: 1D advection operators in FV3. This is designed to be an *interactive* notebook for downloading and playing with the options, initial conditions, zooms, and so on.

fv3_level_transmogrifier
: An *interactive* notebook that shows different hybrid-level setups within FV3, and allows detection of discontinuities within the levels that may cause errors or instabilities. The directory contains the notebook and its dependencies, including source for a Python-wrapped version of `set_eta()`. 

## 2D Global Shallow-water Cases
RHwave
: Rossby-Haurwitz wave, a test of height-vorticity consistency

BTwave
: Barotropic instability, demonstrating vorticity preservation and wave breaking (cf. Galewsky et al. 2004; Scott, Harris, and Polvani 2016)

BLvortex
: Bates-Li forced polar vortex

SWmodon
: Lin-Chen-Yao modon demonstrating the crucial nature of nonlinear vorticity dynamics. Notebook forthcoming; see [Lin et al. (JAMES, 2017)](http://dx.doi.org/10.1002/2017MS000965)

## 3D Global Cases

BCwave
: Hydrostatic baroclinic wave, with and without moisture

TC
: Reed-Jablonowski TC tests, demonstrating the effect of advection schemes and numerical diffusion

mtn_rest_100km
: Resting atmosphere over oscilliatory topography, to diagnose pressure-gradient force truncation error on large scales

TornadicSupercell
: Global super-stretched grid with Toy semi-circle hodograph creating supercell thunderstorms with tornado-like vortices, demonstrating the importance of vorticity preservation at kilometer scales. See animations on [Google Drive](https://drive.google.com/drive/folders/1pVNAuKrYKwxVAlCdVa5faIVRBaK2hdVI) (noaa.gov only).

HSzuritasuperrotation
: The classic Held-Suarez idealized climate, with modifications from [Zurita-Gotor et al. (JAS, 2022)](https://journals.ametsoc.org/view/journals/atsc/79/5/JAS-D-21-0269.1.xml) to demonstrate superrotation. 

## 2D Periodic

MountainWaveIC
: A demonstration of how to rigorously compute thermodynamic quantities in FV3 of importance for mountain wave simulation

mtn_wave_tests
: Standard mountain wave over Schar topography, to demonstrate mountain-wave propagation and 2D FV3 capabilities; and a resting atmosphere over Schar topography, to diagnose pressure-gradient force truncation error on small scales and errors due to hybridization of the vertical coordinate

## 3D Doubly-periodic Cases

DPsupercell
: Supercell on a doubly-periodic using Weisman (WK82) sounding and a straight-line hodograph. Demonstrates pressure partitioning between moist, dry, and nonhydrostatic contributions.
