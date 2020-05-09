ARG NCEP_IMAGE=fv3atm-ncep
ARG FMS_IMAGE=fv3atm-fms
ARG ESMF_IMAGE=fv3atm-esmf

## ---------------------------------------------------------------------------------
## Build environment image
FROM ubuntu:19.10 AS fv3atm-environment

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
FROM fv3atm-environment AS fv3atm-ncep

# download and install NCEP libraries
RUN git config --global http.sslverify false && \
    git clone https://github.com/NCAR/NCEPlibs.git /NCEPlibs && \
    mkdir /opt/NCEPlibs && \
    cd NCEPlibs && \
    git checkout 3da51e139d5cd731c9fc27f39d88cb4e1328212b && \
    echo "y" | ./make_ncep_libs.sh -s linux -c gnu -d /opt/NCEPlibs -o 1

## ---------------------------------------------------------------------------------
## Build ESMF
FROM fv3atm-environment AS fv3atm-esmf

ENV ESMF_DIR=/esmf
ENV ESMF_INSTALL_PREFIX=/usr/local/esmf
ENV ESMF_NETCDF_INCLUDE=/usr/include
ENV ESMF_NETCDF_LIBS="-lnetcdf -lnetcdff"
ENV ESMF_BOPT=O3

RUN git clone -b ESMF_8_0_0 --depth 1 https://git.code.sf.net/p/esmf/esmf $ESMF_DIR && \
    cd $ESMF_DIR && \
    make lib -j8 && \
    make install && \
    make installcheck

## ---------------------------------------------------------------------------------
## Build FMS
FROM fv3atm-environment AS fv3atm-fms

ENV CC=mpicc
ENV FC=mpifort
ENV LDFLAGS="-L/usr/lib"
ENV LOG_DRIVER_FLAGS="--comments"
ENV CPPFLAGS="-I/usr/include -Duse_LARGEFILE -DMAXFIELDMETHODS_=500 -DGFS_PHYS"
ENV FCFLAGS="-fcray-pointer -Waliasing -ffree-line-length-none -fno-range-check -fdefault-real-8 -fdefault-double-8 -fopenmp"

COPY patches/PATCH-FMS-affinity.c-fix-getpid.patch /
COPY patches/PATCH-FMS-libFMS_Makefile.am-add-sat-vapor-pres.patch /
RUN git clone https://github.com/NOAA-GFDL/FMS.git /FMS && \
    cd /FMS && git checkout f68878549364e39591216874d5a922c06b7d3a91 && \
    git am /PATCH-FMS-affinity.c-fix-getpid.patch && \
    git am /PATCH-FMS-libFMS_Makefile.am-add-sat-vapor-pres.patch && \
    autoreconf --install && ./configure && \
    make -j8 && \
    mv /FMS/*/*.mod /FMS/*/*.o /FMS/*/*.h /FMS/

## ---------------------------------------------------------------------------------
## Build fv3atm executable in its own image
FROM $NCEP_IMAGE AS ncep_image
FROM $FMS_IMAGE AS fms_image
FROM $ESMF_IMAGE AS esmf_image
FROM fv3atm-environment AS fv3atm-build

ENV NCEP_DIR=/NCEPlibs \
    FMS_DIR=/FMS \
    ESMF_DIR=/usr/local/esmf

ENV ESMF_INC="-I/usr/local/esmf/include -I${ESMF_DIR}/mod/modO3/Linux.gfortran.64.mpiuni.default/" \
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${ESMF_DIR}/lib/libO3/Linux.gfortran.64.mpiuni.default/:${FMS_DIR}/libFMS/.libs/

COPY --from=ncep_image /opt/NCEPlibs $NCEP_DIR
COPY --from=fms_image /FMS $FMS_DIR
COPY --from=esmf_image /usr/local/esmf $ESMF_DIR

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
FROM fv3atm-environment AS fv3atm-deploy

RUN apt-get update && apt-get install -y \
    python3 \
    python3-pip

RUN ln -s /bin/python3 /bin/python && \
    ln -s /bin/pip3 /bin/pip && \
    pip install --no-cache-dir pyyaml

ENV LD_LIBRARY_PATH=/usr/local/lib
COPY --from=fv3atm-fms /FMS/libFMS/.libs/* /usr/local/lib/
COPY --from=fv3atm-esmf /usr/local/esmf/lib/libO3/Linux.gfortran.64.mpiuni.default/* /usr/local/lib/

RUN git clone -b v0.3.2 https://github.com/VulcanClimateModeling/fv3config.git && \
    cd fv3config && \
    pip install --no-cache-dir -e .

RUN mkdir /work
WORKDIR /work
COPY --from=fv3atm-build /FV3/fv3.exe /FV3/fv3.exe
#COPY etc/create_rundir.py /work
#COPY etc/default.yaml /work
#COPY etc/run_test.sh /work

# run model
CMD ["/work/run_test.sh"]