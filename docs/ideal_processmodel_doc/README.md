Some notes on the ideal "Process Model" setup: 
 
The idealized "Process Model" provides the option to initialize FV3 with horizontally homogeneous environment in a doubly-periodic domain with any physics suite. The environmental setup is done through the fv_processmodel_nml namelist. This option is separate and distinct from the test cases in test_case.F90 that provide specific dycore tests. Here we provide the ability for the user to choose their set up at run time. Options include several different hard-coded thermal and wind profile soundings along with the ability to read in a provided sounding file. Users can also add initial thermal bubbles at their chosen locations and intensities.

Use of this capability require sfc_data.nc and oro_data.nc even if surface/PBL physics are not utilized. The grids in these files must be at least 1 grid point larger, in both x & y directions, than npx & npy. A configurable python script (create_ideal_sfc_oro_input.py) is provided here to generate idealized versions of these files. 

See driver/UFS/fv_processmodel.F90 in the atmos_cubed_sphere directory for code specifics.

A self-contained experiment directory (working as of initial release) can be downloaded from 
  https://github.com/MicroTed/ufs_ideal_run.git

Example namelist settings:

  &fv_processmodel_nml
      bubble_type = 3
	! 0 (default) : no bubble ; 1 : single bubble in domain center ; 2 : Line of N-S bubbles on left side of domain to initialize a squall line ; 3 : User-entered x,y locations of n_nun bubble centers
      n_bub = 3
	! number of bubbles (valid only if bubble_type = 3)
      jcenters = 60, 60, 60
	! j index of bubble centers; must be n_bub in length (valid only if bubble_type = 3)
      icenters = 25, 50, 75
	! i index of bubble centers; must be n_bub in length (valid only if bubble_type = 3)
      t_profile = -1
	! 0 prescribed Weisman-Klemp sounding (default)
	! 1 adiabatic sounding with adi_t K value
	! 2 isothermal sounding with iso_t K value
	! -1 read in sounding from external file input_sounding
      q_profile = -1
	! 0 prescribed Weisman-Klemp sounding (default)
	! 1 dry sounding
	! -1 read in sounding from external file input_sounding
      ws_profile = -1
	! 0 quarter circle hodograph (Harris Approximation) scaled by us0 (default)
	! 1 unidirectional Weisman-Klemp shear profile scaled by us0
	! 2 constant u = us0, v = 0
	! 3 constant v=us0, u=0
	! 4 no wind
	! -1 read in sounding from external file input_sounding
      bubble_t = 2.0 ! max center T (K) perturbation of bubbles
      bubble_q = 0. ! max center qv (g/kg) perturbation 
      bubble_rad_x = 4000. ! bubble radius (m) in the x-direction (m; 10000 default)
      bubble_rad_y = 4000. ! bubble radius (m) in the x-direction (m; 10000 default)
      bubble_zc = 1500. ! bubble radius (m) in the z-direction (m; 1400 default)
      do_coriolis = .false. ! set true for f-plane based using value of deglat
      us0 = 12. ! scaling for wind profiles (m/s ; default 30)
      adi_th = 300. ! adiabatic temperature (K) for t_profile=1 (default 300 K)
      iso_th = 300. ! isothermal temperature (K) for t_profile=2 (default 300 K)
      do_rand_perts = .false. ! if n_bub > 0, applies small amplitude(0.2k; 1E-7 kg/kg), random perturbations to T and qv values inside bubbles (default=.false.) NOTE: These random perts are not currently reproducible under different domain decompositions (i.e., different number of MPI threads)
      umove = 0.0 ! 'grid motion' u-component (m/s) subtracted from sounding U-wind to help keep storm centered in the domain
      vmove = 0.0 ! 'grid motion' v-component (m/s) subtracted from sounding V-wind to help keep storm centered in the domain
  /

Other settings:

  &fv_core_nml
      external_eta = .false. !required .false.
      external_ic = .false. !required .false.
      grid_type = 4      ! selects Cartesian periodic grid; required 4
      npz_type = 'meso' ! options here are buried in the code; ‘meso’ provides a lower top than       !operational configuration (~20 hPa)
      dx_const = 3000. ! set to desired value of dx (meters)
      dy_const = 3000. ! set to desired value of dy (meters)
      regional = .false. ! Required .false.
      deglat = 15 ! grid latitude (default 15 degrees) Affects radiation and coriolis (f-plane)
      deglon = 0  ! grid longitude (default 0) Only affects radiation (time of day relative to UTC)
  /

suite FV3_mp_nssl_ideal:
Simplest physics suite is microphysics only, with this one set up to run the NSSL MP scheme. 

example input.nml: input.nml.idealbubble.test.3m (NSSL 3-moment mp)

suite FV3_ideal_pbl_mp_nssl:

This suite adds PBL/Surface physics. It needs input files
 aerosol.dat
 sfc_emissivity_idx.txt
 solarconstant_noaa_an.txt
 co2historicaldata_*

 Example input.nml: input.nml.idealpblbubble.test.3m (NSSL 3-moment mp)



