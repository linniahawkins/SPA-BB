#!/usr/bin/env python2.7

import numpy as np
import pandas as pd
import os
from scipy import optimize
from pyfunc import *
from datetime import datetime
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters

###############################################################################

def readin_met(in_dir,start,end):
	filename = os.path.join(in_dir+'FLX_USMe2_2002-01-01_2014-12-31_30min.csv')
	df = pd.read_csv(filename,header=0,index_col=0, parse_dates=True, squeeze=True)
	timestamp_hh = pd.date_range(datetime(2002,1,1,0,0,0), datetime(2014,12,31,23,59,59), freq='30min') #'30min','H'
	df.index = timestamp_hh
	met_data = df[start:end]

	return met_data

###############################################################################

def CiFunc(ci_val):

	global gs, et, lt, vcmax, vjmax, kc ,ko, gamma1, resp

	# Photosynthesis
	anx = farquhar ( vcmax, vjmax, kc, ko, gamma1, resp, par, ci_val)

	# calculate co2 at leaf surface
	#cs = co2 - (anx / gbc) # gbc boundary layer conductance to CO2
	cs = co2

	if ( model == 0 ):

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
	elif ( model == 1 ):
		
		# Medlyn
		vpd_min = 1 # (kPa)
		g0 = 0.01
		g1 = 2.8

		if (anx > 0):
			vpd_term = max(vpd,vpd_min) # units kPa
			term = 1.6 * anx / cs
			aquad = 1
			bquad = -(2 * (g0 + term) + (g1 * term)**2 / (gbv * vpd_term))
			cquad = g0 * g0 + (2*g0 + term * (1-g1*g1 / vpd_term)) * term
			[r1,r2] = quadratic(aquad,bquad,cquad)
			gsx = max(r1,r2)
		else:
			gsx = g0

	# convert gs (mol/m2/s) to m/s
	gs = gsx*(Rcon * ( 273.15 + temp )) / atmos_press

	lt = leaf_temperature( gs, netrad, temp, wdef, gbb )

	et = evap( gs, lt, netrad, wdef, gbb )

	adx = diffusion( gs, ci_val, et , gbb , lt, atmos_press, co2)

	diff = adx - max(0, anx)

	return diff

###############################################################################

def empirical_stomata( tleaf, model ):
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

	# modify vcmax and stomatal parameters by Btran
	btran = 1
	vcmax = vcmax * btran

	# calculate leaf respiration at leaf temperature (umol/m2/s)
	resp = rn * nit * np.exp( np.log(2) * ( lt - 10 ) / 10 )

	# calculate saturation vapor pressure and deficit at leaf temperature
	# in progress

	# Ci calculation
	ci0 = 0.1*co2 # minimum plausible Ci
	ci1 = 2*co2 # maximum plausible ci
	tol = .1 # convergence tolerance for brentq solver

	ci = optimize.brentq(CiFunc,ci0,ci1)

	# calculate relative humidity and VPD at leaf surface given new gs
	# in progress

	return lt

###############################################################################

def TleafFunc(tleaf_in):
	""" calculate leaf fluxes from input leaf temperature (tleaf_in) and 
	compare the new temperature to the prior temperature. This function 
	equals zero when tleaf does not change between iterations. """

	global lt

	tleaf_old = tleaf_in

	empirical_stomata(tleaf_old, model)
	tleaf_new = lt

	TleafFunc = tleaf_new - tleaf_old

	return TleafFunc

########################################################
#-------------------- Main code -----------------------#
########################################################

# set data directory
in_dir='/Users/linniahawkins/Documents/SPA/inputs/'

# set simulations start:end dates
start=datetime(2012,7,1) 
end=datetime(2012,7,7)

# set which model to run: 0=Ball-Berry 1=Medlyn
model = 1

# readin met data
met_data = readin_met(in_dir,start,end) 

# set veg parameters
LAI = 3.4 # leaf area index
nit = 2.1 # canopy layer specific nitrogen content (gN/m2 leaf area)
kappac = 16.86
kappaj = 45.37
vcm = kappac * nit # calculate Vcmax/Jmax given current nitrogen (umolC.m-2.s-1)
vjm = kappaj * nit 
dimen = 0.002 # leaf characteristic dimension (m)

# set miscelanious constants
tower_ht = 33 # (m)
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


# initialize output data
out_data = []

# loop through timesteps
for i in range(len(met_data)):

	# set met data for time step
	temp = met_data['airt'][i]
	par = met_data['PAR'][i]
	atmos_press =  98400 # Pa
	netrad =  .0001*met_data['sw_rad'][i] # net radiation penetration to canopy layer (kW.m-2)
	co2 = met_data['co2'][i] # umol/mol
	vpd = met_data['vpd'][i] # kPa
	wind_spd = met_data['wind_spd'][i]  # m/s

	# Calculate relative humidity (%) and water deficit (g/m3) from VPD (kPa) and temp (C)
	if (vpd<0):
		vpd = 0
	wdef = vpd * 217 / (.1 * (temp + 273.4))  # air water content deficit (g m-3)
	vpsat = 613.75*np.exp((17.502*temp) / (240.97+temp)) # (Pa)
	vpair = vpsat - vpd*(10**3) # (Pa)
	# constrain vpair
	vpair = max(min(vpair,0.2*vpsat),vpsat)
	rh = (vpair/vpsat)

	[gbb, gbh] = boundary(temp,tower_ht, atmos_press, wind_spd, dimen )
	gbv = gbb*atmos_press/(Rcon*(temp+273.15)) # boundary layer conductance to H2O (umol/m2/s)


	# Run Stomatal Optimization Model 
	lt_out = optimize.brentq(TleafFunc, temp-20, temp+20)
	

	# prepare output data
	an = farquhar ( vcmax, vjmax, kc, ko, gamma1, resp, par, ci)
	agr = an+resp
	gsm = 1000*gs* atmos_press / ( ( lt + 273.15 ) * Rcon ) # convert m/s to mmol/m2/s
	et = et/1000 # convert from g/m2/s to kg/m2/s
	
	# for comparison with SPA-WUEi scale by leaf area: la = .16667 * 3.4 = 1/6 * LAI 
	la = LAI/6
	an_total = an * la
	agr_total = agr * la
	resp_total = resp * la
	et_total = et * la 
	gs_total = gsm * la

	# store data for this time step
	out_data.append(
		{
			'gs(mmol/m2/s)'	:	gs_total,
			'tleaf(C)'		:	lt,
			'anet(umol/m2/s)':	an_total,
			'agr(umol/m2/s)': 	agr_total,
			'resp(umol/m2/s)': 	resp_total,
			'ci(umol/mol)'	:	ci,
			'et(kg/m2/s)'	:	et_total,

		})

# convert out_data to dataframe
out_data = pd.DataFrame(out_data)
out_data.index = pd.date_range(start,end,freq='30min')

# Write output data
#out_file = '/Users/linniahawkins/Documents/SPA/spa-bb/python_bb/SPA_Med_python_1layer_' + str(datetime.date(start)) + '_' + str(datetime.date(end)) + '.csv' 
#out_data.to_csv(out_file)

########################################################
#---------------------- Plot --------------------------#
########################################################
plt.figure(num=None, figsize=(12, 9), dpi=100, facecolor='w', edgecolor='k')
plt.subplots_adjust(left=.05, bottom=.1, right=.9, top=.95, wspace=0.2, hspace=0.2)

ax1 = plt.subplot(3,1,1)
ax1.plot(out_data['anet(umol/m2/s)'],label='Anet (umol/m2/s)')
ax1.set_ylabel('(umol/m2/s)')
ax2 = ax1.twinx()
ax2.plot(out_data['et(kg/m2/s)']*86400,color='grey',alpha=0.8,label='ET (mm/day)')
ax2.set_ylabel('(mm/day)')
ax2.legend(loc='upper right',frameon=False)
ax1.legend(loc='upper left',frameon=False)

ax1 = plt.subplot(3,1,2)
ax1.plot(out_data['ci(umol/mol)'],label='Ci (umol/mol)')
ax1.plot(met_data['co2'],label='CO2(umol/mol)')
ax1.set_ylabel('(umol/mol)')
ax1.set_ylim([0,450])
ax2 = ax1.twinx()
ax2.plot(out_data['gs(mmol/m2/s)'],color='grey',alpha=0.8,label='gs (mmol/m2/s)')
ax2.set_ylabel('(mmol/m2/s)')
ax2.legend(loc='upper right',frameon=False)
ax1.legend(loc='upper left',frameon=False)

ax1 = plt.subplot(3,1,3)
ax1.plot(out_data['tleaf(C)'],label='Tleaf (C)')
ax1.plot(met_data['airt'], label='Tair (C)')
ax1.set_ylabel('(C)')
ax2 = ax1.twinx()
ax2.plot(met_data['vpd'],color='grey',alpha=0.8,label='VPD (hPa)')
ax2.set_ylabel('(hPa)')
ax2.legend(loc='upper right',frameon=False)
ax1.legend(loc='upper left',frameon=False)


out_figure = '/Users/linniahawkins/Documents/SPA/spa-bb/python_bb/SPA_MED_python_' + str(datetime.date(start)) + '_' + str(datetime.date(end))
plt.savefig(out_figure, dpi=300)
plt.show()

