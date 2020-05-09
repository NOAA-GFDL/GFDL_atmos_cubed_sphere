#!/usr/bin/env python

import fv3config

fv3config.ensure_data_is_downloaded()
config = fv3config.config_from_yaml('default.yaml')
fv3config.write_run_directory(config, './rundir')

