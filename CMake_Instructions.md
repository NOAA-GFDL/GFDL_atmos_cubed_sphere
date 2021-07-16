Load the compiler, mpi, netCDF and FMS module

```sh
git clone -b feature/cmake_in_dycore_of_emc https://github.com/noaa-emc/GFDL_atmos_cubed_sphere fv3dycore_emc
mkdir build && cd build
cmake -DGFS_PHYS=ON ../fv3dycore_emc
make -j8
```

The `cmake` configuration will succeed, but the `make` will fail due to unresolved CCPP dependencies.
It is designed to fail due to CCPP dependencies.
