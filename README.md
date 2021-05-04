# GFDL_atmos_cubed_sphere

This is for the FV3 dynamical core and the GFDL Microphysics for use by NCEP/EMC within GFS.
The source in this branch reflects the codebase used by NCEP/EMC for use in GFS and UFS.
# Where to find information

Visit the [FV3 website](https://www.gfdl.noaa.gov/fv3/) for more information. Reference material is available at [FV3 documentation and references](https://www.gfdl.noaa.gov/fv3/fv3-documentation-and-references/). 

# Proper usage attribution

Cite Putman and Lin (2007) and Harris and Lin (2013) when describing a model using the FV3 dynamical core.
Cite Chen et al (2013) and Zhou et al (2019) when using the GFDL Microphysics.

# What files are what

The top level directory structure groups source code and input files as follow:

| File/directory       | Purpose |
| --------------       | ------- |
| ```LICENSE.md```     | a copy of the Gnu lesser general public license, version 3. |
| ```README.md```      | this file with basic pointers to more information |
| ```model/```         | contains the source code for core of the FV3 dyanmical core |
| ```driver/```        | contains drivers used by different models/modeling systems |
| ```tools/```         | contains source code of tools used within the core |
| ```docs/```          | contains documentation for the FV3 dynamical core |

# Disclaimer

The United States Department of Commerce (DOC) GitHub project code is provided
on an "as is" basis and the user assumes responsibility for its use. DOC has
relinquished control of the information and no longer has responsibility to
protect the integrity, confidentiality, or availability of the information. Any
claims against the Department of Commerce stemming from the use of its GitHub
project will be governed by all applicable Federal law. Any reference to
specific commercial products, processes, or services by service mark,
trademark, manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce. The
Department of Commerce seal and logo, or the seal and logo of a DOC bureau,
shall not be used in any manner to imply endorsement of any commercial product
or activity by DOC or the United States Government.

This project code is made available through GitHub but is managed by NOAA-GFDL
at https://gitlab.gfdl.noaa.gov.
