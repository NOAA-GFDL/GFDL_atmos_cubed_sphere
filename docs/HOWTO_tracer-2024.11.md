NOTE: these tracers are not specific ratios and so should not be mass-adjusted.  
To activate, add the following to your field_table:  

```
 "TRACER", "atmos_mod", "w_diff"
           "longname",     "w_diff"
           "units",        "m/s"
           "adjust_mass", "false"
       "profile_type", "fixed", "surface_value=0" /
 "TRACER", "atmos_mod", "pbl_age"
           "longname",    "Age of air from PBL"
           "units",     "d"
           "adjust_mass", "false"
       "profile_type", "fixed", "surface_value=0." /
 "TRACER", "atmos_mod", "tro_pbl_age"
           "longname",    "Age of air from tropical PBL"
           "units",     "d"
           "adjust_mass", "false"
       "profile_type", "fixed", "surface_value=0." /
```
