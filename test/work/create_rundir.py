#!/usr/bin/env python

import fv3config

rundir = './rundir'
config = './default.yaml'

print('Creating rundir [' + rundir + '] from config [' + config + ']')
config = fv3config.config_from_yaml(config)
fv3config.write_run_directory(config, rundir)

