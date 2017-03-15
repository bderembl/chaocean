#!/usr/bin/env python

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt
import pupynere as netcdf
import MITgcmutils as mit

binprec = '>f4'


dir0 = '/tank/chaocean/run05/mnc*/'
file1 = 'state.0000*'

dir_o = '/tank/groups/climode/chaocean/restart/';

f1 = mit.mnc_files(dir0 + file1)

T = f1.variables['T'][:]
nt = len(T)-1

uvel  = f1.variables['U'   ][nt,:,:,:]
vvel  = f1.variables['V'   ][nt,:,:,:]
theta = f1.variables['Temp'][nt,:,:,:]
salt  = f1.variables['S'   ][nt,:,:,:]
eta   = f1.variables['Eta' ][nt,:,:]


si_z,si_y,si_x = theta.shape

uvel = uvel[:si_z,:si_y,:si_x]
vvel = vvel[:si_z,:si_y,:si_x]

uvel.astype(binprec).tofile( dir_o + 'uinit.box')
vvel.astype(binprec).tofile( dir_o + 'vinit.box')
theta.astype(binprec).tofile(dir_o + 'tinit.box')
salt.astype(binprec).tofile( dir_o + 'sinit.box')
eta.astype(binprec).tofile(  dir_o + 'einit.box')
