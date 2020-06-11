#!/usr/bin/env python2.7

import netCDF4 as nc
import numpy as np
import pandas as pd
import os
from scipy import optimize
from pyfunc import *

# Set atmospheric variables 
rh = 65 # %
temp = 25 #C
par = 790
atmos_press =  98400 # Pa
lt = temp +1
netrad =  .02 # net radiation penetration to canopy layer (kW.m-2)
co2 = 385 # umol/mol

# set miscelanious values
Rcon = 8.3144    # constant
gbb = 0.002     # boundary layer conductance m/s
gbv = gbb * (10**6) / 18 # boundary layer conductance H2O (mol H2O/m2 leaf/s)
rn = 0.105 # respiration constant umol CO2/g N at 10 deg C
nit = 2.5 # canopy layer specific nitrogen content (gN/m2 leaf area)


# Set photosynthesis parameters
vcmax = 18
vjmax = 30
kc = 40
ko = 45
gamma1 = 16

wdef = 10
lt = temp
resp = rn * nit * np.exp( np.log(2) * ( lt - 10 ) / 10 )

def CiFunc(ci_val):

    # Photosynthesis
    anx = farquhar ( vcmax, vjmax, kc, ko, gamma1, resp, par, ci_val)

    # calculate co2 at leaf surface
    cs = co2 - (anx / gbv)
   
    # Ball-Berry
    g0 = 0.01
    g1 = 9

    aquad = cs
    bquad = cs*(gbv-g0 + g1*anx*rh )
    cquad = -gbv * (cs*g0 + g1*anx*rh)
    [r1,r2] = quadratic(aquad,bquad,cquad)
    gsx = max(r1,r2)

    # convert gs (mol/m2/s) to m/s
    convert = (Rcon * ( 273.15 + temp )) / atmos_press
    gs = gsx*convert
    #gs = 0.002 # m/s

    lt = leaf_temperature( gs, netrad, temp, wdef, gbb )

    et = evap( gs, lt, netrad, wdef, gbb )
    #et = 0.001

    adx = diffusion( gs, ci_val, et , gbb , lt, atmos_press, co2)

    diff = adx - max(0, anx)

    return diff


##### Main code #####

ci0 = 0.5*co2
ci1 = 2*co2
tol = 0.1
ci = optimize.brent(CiFunc,brack=(ci0,ci1),tol=tol)

##### plot #####
print ci

