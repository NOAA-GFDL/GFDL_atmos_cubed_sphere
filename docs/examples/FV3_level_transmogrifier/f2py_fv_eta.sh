rm -rf *.o *.mod *.so
f2py -c -m fv_eta fv_eta.F90 -I/Users/wolfganglanghans/wavespotter/GFDL_atmos_cubed_sphere/docs/examples/FV3_level_transmogrifier
