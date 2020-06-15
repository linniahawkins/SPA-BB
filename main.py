#!/usr/bin/env python2.7

import numpy as np
import pandas as pd
import os
from scipy import optimize
from pyfunc import *

# Set atmospheric variables 
rh = 65 # %
temp = 25 #C
par = 700
atmos_press =  98400 # Pa
netrad =  .02 # net radiation penetration to canopy layer (kW.m-2)
co2 = 385 # umol/mol

# set veg parameters
la = 3.2
nit = 2.1 # canopy layer specific nitrogen content (gN/m2 leaf area)
kappac = 16.86
kappaj = 45.37

# calculate Vcmax/Jmax given current nitrogen (umolC.m-2.s-1)
vcm = kappac * nit
vjm = kappaj * nit 


# set miscelanious values
Rcon = 8.3144    # constant
rn = 0.105 # respiration constant umol CO2/g N at 10 deg C
Vc_kurtosis = 0.143 # kurtosis of Vcmax temp response
Vj_kurtosis = 0.172 # Kurtosis of Jmax temp response
kc_saturation = 310
kc_half_sat_conc = 23.956
ko_saturation = 155
ko_half_sat_conc = 14.509
co2comp_saturation = 36.5
co2comp_half_sat_conc = 9.46
vcmax_max_temp = 65.03 # max tolerated temperature for carboxylation
jmax_max_temp = 57.05 # max tolerated temperature for electron transport
metabolic_opt_temp = 30 # metabolic temperature optimum (C)



gbb = 0.002     # boundary layer conductance m/s
gbv = gbb * (10**6) / 18 # boundary layer conductance H2O (mol H2O/m2 leaf/s)
#gbv =  1 mol/m2/s 



lt = temp

# Apply temperature modifications to metabolic parameters
vcmt   = tempmet(vcmax_max_temp, metabolic_opt_temp, Vc_kurtosis, lt)
vjmt   = tempmet(jmax_max_temp,metabolic_opt_temp,Vj_kurtosis,lt)
vcmax  = vcmt * vcm
vjmax  = vjmt * vjm
kc     = arrhenious(kc_saturation,kc_half_sat_conc,lt)
ko     = arrhenious(ko_saturation,ko_half_sat_conc,lt)
gamma1 = arrhenious(co2comp_saturation,co2comp_half_sat_conc,lt)


wdef = 10
resp = rn * nit * np.exp( np.log(2) * ( lt - 10 ) / 10 )

def CiFunc(ci_val):

	global gs, lt, et
	print ci_val

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
	print gs

	lt = leaf_temperature( gs, netrad, temp, wdef, gbb )
	print lt

	et = evap( gs, lt, netrad, wdef, gbb )

	adx = diffusion( gs, ci_val, et , gbb , lt, atmos_press, co2)

	diff = adx - max(0, anx)


	return diff


def empirical_stomata( tleaf ):
	"""determine stable ci for given leaf temperature
		tleaf	:	leaf temperature (C)"""

	global ci, lt

	# Apply temperature modifications to metabolic parameters
	vcmt   = tempmet(vcmax_max_temp, metabolic_opt_temp, Vc_kurtosis, tleaf)
	vjmt   = tempmet(jmax_max_temp,metabolic_opt_temp,Vj_kurtosis,tleaf)
	vcmax  = vcmt * vcm
	vjmax  = vjmt * vjm
	kc     = arrhenious(kc_saturation,kc_half_sat_conc,tleaf)
	ko     = arrhenious(ko_saturation,ko_half_sat_conc,tleaf)
	gamma1 = arrhenious(co2comp_saturation,co2comp_half_sat_conc,tleaf)

	# modify vcmax fo Btran

	# calculate leaf respiration at leaf temperature (umol/m2/s)
	resp = rn * nit * np.exp( np.log(2) * ( tleaf - 10 ) / 10 )

	# calculate saturation vapor pressure and deficit at leaf temperature

	# Ci calculation
	ci0 = 0.1*co2 # minimum plausible Ci
	ci1 = co2 # maximum plausible ci
	tol = 0.1 # convergence tolerance for brentq solver

	ci = optimize.brentq(CiFunc,ci0,ci1)

	# calculate relative humidity and VPD at leaf surface given new gs
	#wdef = ?

	return lt



def TleafFunc(tleaf_in):
	""" calculate leaf fluxes from input leaf temperature (tleaf_in) and 
	compare the new temperature to the prior temperature. This function 
	equals zero when tleaf does not change between iterations. """

	global lt 

	delta = 0.05 

	lt = empirical_stomata(tleaf_in-delta)
	tleaf_old = lt
	lt = empirical_stomata(tleaf_in)
	tleaf_new = lt
	TleafFunc = tleaf_new - tleaf_old
	print TleafFunc

	return TleafFunc



##### Main code #####

lt_out = optimize.brentq(TleafFunc,temp-20,temp+20,xtol = 0.1)

print gs, ci, lt


"""
ci0 = 0.1*co2 # minimum plausible Ci
ci1 = co2 # maximum plausible ci
tol = 0.1 # convergence tolerance for brentq solver

ci_out = optimize.brentq(CiFunc,ci0,ci1)
an = farquhar ( vcmax, vjmax, kc, ko, gamma1, resp, par, ci_out)
"""




