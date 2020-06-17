#!/usr/bin/env python2.7

import numpy as np
import pandas as pd
import os
from scipy import optimize
from pyfunc import *
from datetime import datetime
import matplotlib.pyplot as plt

# Set atmospheric variables 
def readin_met(in_dir,start,end):
	filename = os.path.join(in_dir+'FLX_USMe2_2002-01-01_2014-12-31_30min.csv')
	df = pd.read_csv(filename,header=0,index_col=0, parse_dates=True, squeeze=True)
	timestamp_hh = pd.date_range(datetime(2002,1,1,0,0,0), datetime(2014,12,31,23,59,59), freq='30min') #'30min','H'
	df.index = timestamp_hh
	met_data = df[start:end]

	return met_data


in_dir='/Users/linniahawkins/Documents/SPA/inputs/'  
start=datetime(2012,7,1) 
end=datetime(2012,7,31)

met_data = readin_met(in_dir,start,end) 
varb = ['TA_F','SW_IN_F','VPD_F','P_F','SWC_F_MDS_1']

met_temp = met_data['airt'] #C
met_par = met_data['PAR']
atmos_press =  98400 # Pa
met_netrad =  met_data['sw_rad']*.1 # net radiation penetration to canopy layer (kW.m-2)
met_co2 = met_data['co2'] # umol/mol
met_vpd = met_data['vpd'] # kPa
met_wind_spd = met_data['wind_spd'] # m/s


# set veg parameters
la = 3.2
nit = 2.1 # canopy layer specific nitrogen content (gN/m2 leaf area)
kappac = 16.86
kappaj = 45.37
vcm = kappac * nit # calculate Vcmax/Jmax given current nitrogen (umolC.m-2.s-1)
vjm = kappaj * nit 
dimen = 0.002 # (m)
tower_ht = 33 # (m)


# set miscelanious values
Rcon = 8.3144    # universal gas constant (J/K/mol)
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


def CiFunc(ci_val):

	global gs, et, lt, vcmax, vjmax, kc ,ko, gamma1, resp

	# Photosynthesis
	anx = farquhar ( vcmax, vjmax, kc, ko, gamma1, resp, par, ci_val)

	# calculate co2 at leaf surface
	#cs = co2 - (anx / gbc) # gbc boundary layer conductance to CO2
	cs = co2

	# Ball-Berry
	g0 = 0.01
	g1 = 9

	if (anx > 0 ):
		aquad = cs
		bquad = cs*(gbv-g0) - ( g1*anx )
		cquad = -gbv * (cs*g0 + g1*anx*rh)
		[r1,r2] = quadratic(aquad,bquad,cquad)
		gsx = max(r1,r2)
	else:
		gsx = g0

	# convert gs (mol/m2/s) to m/s
	convert = (Rcon * ( 273.15 + temp )) / atmos_press
	gs = gsx*convert

	lt = leaf_temperature( gs, netrad, temp, wdef, gbb )

	et = evap( gs, lt, netrad, wdef, gbb )

	adx = diffusion( gs, ci_val, et , gbb , lt, atmos_press, co2)

	diff = adx - max(0, anx)

	return diff


def empirical_stomata( tleaf ):
	"""determine stable ci for given leaf temperature
		tleaf	:	leaf temperature (C)"""

	global gs, et, lt, ci, vcmax, vjmax, kc ,ko, gamma1, resp, wdef

	lt = tleaf

	# Apply temperature modifications to metabolic parameters
	vcmt   = tempmet(vcmax_max_temp, metabolic_opt_temp, Vc_kurtosis, lt)
	vjmt   = tempmet(jmax_max_temp,metabolic_opt_temp,Vj_kurtosis,lt)
	vcmax  = vcmt * vcm
	vjmax  = vjmt * vjm
	kc     = arrhenious(kc_saturation,kc_half_sat_conc,lt)
	ko     = arrhenious(ko_saturation,ko_half_sat_conc,lt)
	gamma1 = arrhenious(co2comp_saturation,co2comp_half_sat_conc,lt)

	# modify vcmax by Btran

	# calculate leaf respiration at leaf temperature (umol/m2/s)
	resp = rn * nit * np.exp( np.log(2) * ( lt - 10 ) / 10 )

	# calculate saturation vapor pressure and deficit at leaf temperature

	# Ci calculation
	ci0 = 0.1*co2 # minimum plausible Ci
	ci1 = 2*co2 # maximum plausible ci
	tol = .1 # convergence tolerance for brentq solver

	ci = optimize.brentq(CiFunc,ci0,ci1)

	# calculate relative humidity and VPD at leaf surface given new gs

	return lt

def TleafFunc(tleaf_in):
	""" calculate leaf fluxes from input leaf temperature (tleaf_in) and 
	compare the new temperature to the prior temperature. This function 
	equals zero when tleaf does not change between iterations. """

	global lt

	delta = 1

	tleaf_old = tleaf_in

	empirical_stomata(tleaf_in)
	#tleaf_old = lt
	#empirical_stomata(tleaf_in)
	tleaf_new = lt

	TleafFunc = tleaf_new - tleaf_old

	return TleafFunc


##### Main code #####



GS=[]; LT=[]; CI=[]; ET=[]; AN=[];

for i in range(len(met_temp)):
	temp = met_temp[i]
	par = met_par[i]
	atmos_press =  98400 # Pa
	netrad =  .08 # net radiation penetration to canopy layer (kW.m-2)
	co2 = met_co2[i] # umol/mol
	vpd = met_vpd[i] # kPa
	wind_spd = met_wind_spd[i]  # m/s

	# Calculate relative humidity (%) and water deficit (g/m3) from VPD (kPa) and temp (C)
	if (vpd<0):
		vpd = 0
	wdef = vpd * 217 / (.1 * (temp + 273.4))  # air water content deficit (g m-3)
	vpsat = 613.75*np.exp((17.502*temp) / (240.97+temp)) # (Pa)
	vpair = vpsat - vpd*(10**3) # (Pa)
	rh = (vpair/vpsat)

	[gbb, leaf_heat_conductance] = boundary(temp,tower_ht, atmos_press, wind_spd, dimen )
	gbv = gbb*atmos_press/(Rcon*(temp+273.15)) # boundary layer conductance to H2O (umol/m2/s)


	lt_out = optimize.brentq(TleafFunc, temp-20, temp+20)
	an = farquhar ( vcmax, vjmax, kc, ko, gamma1, resp, par, ci)
	GS.append(gs)
	LT.append(lt)
	CI.append(ci)
	ET.append(et/1000*86400) # mm/day
	AN.append(an)

############### plot ###################
# scatter vars of interest
plt.figure(num=None, figsize=(12, 9), dpi=100, facecolor='w', edgecolor='k')
plt.subplots_adjust(left=.05, bottom=.1, right=.9, top=.95, wspace=0.2, hspace=0.2)
plt.subplot(2,1,1)
plt.plot(AN,label='an')
plt.plot(ET,label='et')
plt.legend()

plt.subplot(2,1,2)
plt.plot(CI,label='ci')
plt.plot(LT,label='lt')
plt.legend()

plt.show()

