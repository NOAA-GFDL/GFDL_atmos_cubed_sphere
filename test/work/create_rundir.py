#!/usr/bin/env python

import fv3config

config = fv3config.config_from_yaml('default.yaml')
fv3config.write_run_directory(config, './rundir')

