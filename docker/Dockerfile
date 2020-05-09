## ---------------------------------------------------------------------------------
## Build environment image
FROM ubuntu:19.10 AS base

RUN apt-get update && apt-get install -y \
    curl \
    wget \
    gcc \
    git \
    libblas-dev \
    liblapack-dev \
    libnetcdf-dev \
    libnetcdff-dev \
    perl \
    make \
    rsync \
    libffi-dev \
    openssl \
    libopenmpi3 \
    bats

RUN git config --global user.email "you@example.com" && \
    git config --global user.name "Your Name"

## ---------------------------------------------------------------------------------
## Build NCEPlibs
RUN git config --global http.sslverify false && \
    git clone https://github.com/NCAR/NCEPlibs.git /NCEPlibs && \
    mkdir /opt/NCEPlibs && \
    cd NCEPlibs && \
    git checkout 3da51e139d5cd731c9fc27f39d88cb4e1328212b && \
    echo "y" | ./make_ncep_libs.sh -s linux -c gnu -d /opt/NCEPlibs -o 1

## ---------------------------------------------------------------------------------
## Build ESMF
ENV ESMF_DIR=/esmf \
    ESMF_INSTALL_PREFIX=/usr/local/esmf \
    ESMF_NETCDF_INCLUDE=/usr/include \
    ESMF_NETCDF_LIBS="-lnetcdf -lnetcdff" \
    ESMF_BOPT=O3
RUN git clone -b ESMF_8_0_0 --depth 1 https://git.code.sf.net/p/esmf/esmf $ESMF_DIR && \
    cd $ESMF_DIR && \
    make lib -j8 && \
    make install && \
    make installcheck

## ---------------------------------------------------------------------------------
## Build FMS
COPY patches/PATCH-FMS-affinity.c-fix-getpid.patch /
COPY patches/PATCH-FMS-libFMS_Makefile.am-add-sat-vapor-pres.patch /
RUN git clone https://github.com/NOAA-GFDL/FMS.git /FMS && \
    cd /FMS && git checkout f68878549364e39591216874d5a922c06b7d3a91 && \
    git am /PATCH-FMS-affinity.c-fix-getpid.patch && \
    git am /PATCH-FMS-libFMS_Makefile.am-add-sat-vapor-pres.patch && \
    CC=mpicc FC=mpifort LDFLAGS="-L/usr/lib" LOG_DRIVER_FLAGS="--comments" \
      CPPFLAGS="-I/usr/include -Duse_LARGEFILE -DMAXFIELDMETHODS_=500 -DGFS_PHYS" \
      FCFLAGS="-fcray-pointer -Waliasing -ffree-line-length-none -fno-range-check -fdefault-real-8 -fdefault-double-8 -fopenmp" \
      autoreconf --install && \
    CC=mpicc FC=mpifort LDFLAGS="-L/usr/lib" LOG_DRIVER_FLAGS="--comments" \
      CPPFLAGS="-I/usr/include -Duse_LARGEFILE -DMAXFIELDMETHODS_=500 -DGFS_PHYS" \
      FCFLAGS="-fcray-pointer -Waliasing -ffree-line-length-none -fno-range-check -fdefault-real-8 -fdefault-double-8 -fopenmp" \
      ./configure && \
    CC=mpicc FC=mpifort LDFLAGS="-L/usr/lib" LOG_DRIVER_FLAGS="--comments" \
      CPPFLAGS="-I/usr/include -Duse_LARGEFILE -DMAXFIELDMETHODS_=500 -DGFS_PHYS" \
      FCFLAGS="-fcray-pointer -Waliasing -ffree-line-length-none -fno-range-check -fdefault-real-8 -fdefault-double-8 -fopenmp" \
      make -j8 && \
    mv /FMS/*/*.mod /FMS/*/*.o /FMS/*/*.h /FMS/

## ---------------------------------------------------------------------------------
## Build fv3atm executable in its own image
FROM base AS build

ENV NCEP_DIR=/NCEPlibs \
    FMS_DIR=/FMS \
    ESMF_DIR=/usr/local/esmf

ENV ESMF_INC="-I/usr/local/esmf/include -I${ESMF_DIR}/mod/modO3/Linux.gfortran.64.mpiuni.default/" \
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ESMF_DIR}/lib/libO3/Linux.gfortran.64.mpiuni.default/:${FMS_DIR}/libFMS/.libs/

COPY --from=base /opt/NCEPlibs $NCEP_DIR
COPY --from=base /FMS $FMS_DIR
COPY --from=base /usr/local/esmf $ESMF_DIR

ARG configure_file=configure.fv3.gnu_docker

COPY patches/PATCH-stochastic_physics-update_ca.F90-fix-logical-comparisons.patch /
RUN git clone https://github.com/noaa-psd/stochastic_physics.git /stochastic_physics && \
    cd /stochastic_physics && git checkout 01739da377d603c9c8d6d666afe2442260fd79d1

COPY patches/PATCH-fv3atm-apply-fixes.patch /
RUN git clone https://github.com/NOAA-EMC/fv3atm.git /FV3 && \
    cd /FV3 && git checkout b48b640f4d8cae978e2d98e149412054d74dcce7 && \
    git submodule update --init && \
    git am /PATCH-fv3atm-apply-fixes.patch && \
    mkdir namphysics && mkdir conf

RUN /bin/rm -f /PATCH-*
COPY etc/make.rules /FV3/conf/make.rules
COPY etc/configure.fv3.gnu_docker /FV3/conf/configure.fv3

RUN cd /FV3 && make clean && make libs -j8
RUN cd /FV3/atmos_cubed_sphere && make clean && cd /FV3 && make -j8

## ---------------------------------------------------------------------------------
## Build deployment image
FROM base AS test

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip

RUN ln -s /bin/python3 /bin/python && \
    ln -s /bin/pip3 /bin/pip && \
    pip install --no-cache-dir pyyaml

ENV LD_LIBRARY_PATH=/usr/local/lib
COPY --from=fv3atm-base /FMS/libFMS/.libs/* /usr/local/lib/
COPY --from=fv3atm-base /usr/local/esmf/lib/libO3/Linux.gfortran.64.mpiuni.default/* /usr/local/lib/

RUN git clone -b v0.3.2 https://github.com/VulcanClimateModeling/fv3config.git && \
    cd fv3config && \
    pip install --no-cache-dir -e .

RUN mkdir /work
WORKDIR /work
COPY --from=fv3atm-build /FV3/fv3.exe /FV3/fv3.exe

# run model
CMD ["/work/run_test.sh"]