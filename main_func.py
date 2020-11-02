#!/usr/bin/env python2.7

import numpy as np
import pandas as pd
import os
from scipy import optimize
from pyfunc import *
from dics import *
from datetime import datetime
import matplotlib.pyplot as plt
from pandas.plotting import register_matplotlib_converters

###############################################################################

def readin_met(in_file,start,end):
	df = pd.read_csv(in_file,header=0,index_col=0, parse_dates=True, squeeze=True)
	timestamp_hh = pd.date_range(datetime(2002,1,1,0,0,0), datetime(2014,12,31,23,59,59), freq='30min') #'30min','H'
	df.index = timestamp_hh
	met_data = df[start:end]

	return met_data

def met_timestep(met_data):
		# set met data for time step
		temp = met_data['airt']
		par = met_data['PAR']
		atmos_press =  98400 # Pa
		netrad =  .0001*met_data['sw_rad'] # net radiation penetration to canopy layer (kW.m-2)
		co2 = met_data['co2'] # umol/mol
		vpd = met_data['vpd'] # kPa
		wind_spd = met_data['wind_spd']  # m/s

		# Calculate relative humidity (%) and water deficit (g/m3) from VPD (kPa) and temp (C)
		if (vpd<0):
			vpd = 0
		wdef = vpd * 217 / (.1 * (temp + 273.4))  # air water content deficit (g m-3)
		vpsat = 613.75*np.exp((17.502*temp) / (240.97+temp)) # (Pa)
		vpair = vpsat - vpd*(10**3) # (Pa)
		vpair = min(max(vpair,0.05*vpsat),vpsat) # constrain vpair to prevent solution from blowing up
		rh = (vpair/vpsat)

		[gbb, gbh] = boundary(temp,tower_ht, atmos_press, wind_spd, dimen )
		gbv = gbb*atmos_press/(Rcon*(temp+273.15)) # boundary layer conductance to H2O (umol/m2/s)

		met_timestep = {'temp' : temp,
					'par' : par,
					'atmos_press' : atmos_press,
					'netrad' : netrad,
					'co2' : co2,
					'vpd' : vpd,
					'wind_spd' : wind_spd,
					'wdef' : wdef,
					'vpair' : vpair,
					'vpsat' : vpsat,
					'rh' : rh,
					'gbb' : gbb,
					'gbv' : gbv,
					'vpair' : vpair,
					'vpsat' : vpsat,
					}

		return met_timestep


###############################################################################

def CiFunc(ci_val):

	global gs, et, lt, vcmax, vjmax, kc ,ko, gamma1, resp, met_ts

	# modify vcmax and stomatal parameters by Btran
	btran = 1
	vcmax = vcmax * btran

	# Photosynthesis
	anx = farquhar ( vcmax, vjmax, kc, ko, gamma1, resp, met_ts['par'], ci_val)

	# calculate co2 at leaf surface
	#cs = co2 - (anx / gbc) # gbc boundary layer conductance to CO2
	cs = met_ts['co2']


	if ( mod == 0 ):

		# Ball-Berry
		g0 = 0.01
		g1 = g1_param

		if (anx > 0 ):
			aquad = cs
			bquad = cs*(met_ts['gbv']-g0) - ( g1*anx )
			cquad = -met_ts['gbv'] * (cs*g0 + g1*anx*met_ts['rh'])
			[r1,r2] = quadratic(aquad,bquad,cquad)
			gsx = max(r1,r2)
		else:
			gsx = g0

	elif ( mod == 1 ):
		
		# Medlyn
		vpd_min = .1 # (kPa)
		g0 = 0.01
		g1 = g1_param

		if (anx > 0):
			vpd_term = max(met_ts['vpd'],vpd_min) # units kPa
			term = 1.6 * anx / cs
			aquad = 1
			bquad = -(2 * (g0 + term) + (g1 * term)**2 / (met_ts['gbv'] * vpd_term))
			cquad = g0 * g0 + (2*g0 + term * (1-g1*g1 / vpd_term)) * term
			[r1,r2] = quadratic(aquad,bquad,cquad)
			gsx = max(r1,r2)
		else:
			gsx = g0

	# convert gs (mol/m2/s) to m/s
	gs = gsx*(Rcon * ( 273.15 + met_ts['temp'] )) / met_ts['atmos_press']

	lt = leaf_temperature( gs, met_ts['netrad'], met_ts['temp'], met_ts['wdef'], met_ts['gbb'] )

	et = evap( gs, lt, met_ts['netrad'], met_ts['wdef'], met_ts['gbb'] )
	et = max(0,et) / 18 # convert grams/m2/s to mol/m2/s

	adx = diffusion( gs, ci_val, et , met_ts['gbb'] , lt, met_ts['atmos_press'], met_ts['co2'])

	diff = adx - max(0, anx)

	return diff

###############################################################################

def empirical_stomata( tleaf, mod , met_ts):
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


	# calculate leaf respiration at leaf temperature (umol/m2/s)
	resp = rn * nit * np.exp( np.log(2) * ( lt - 10 ) / 10 )

	# calculate saturation vapor pressure and deficit at leaf temperature
	# in progress

	# Ci calculation
	ci0 = 0.1*met_ts['co2'] # minimum plausible Ci
	ci1 = 2*met_ts['co2'] # maximum plausible ci
	tol = .1 # convergence tolerance for brentq solver

	ci = optimize.brentq(CiFunc,ci0,ci1)

	# calculate relative humidity and VPD at leaf surface given new gs
	#  rh_leaf = (gbv*met_ts['vpair'] + gs*met_ts['vpsat']) / ((gbv+gs*met_ts['vpsat']) # fraction
    #  vpd_leaf = max(met_ts['vpsat'] - rh_leaf*met_ts['vpsat'], 0.1) # Pa


	return lt

###############################################################################

def TleafFunc(tleaf_in):
	""" calculate leaf fluxes from input leaf temperature (tleaf_in) and 
	compare the new temperature to the prior temperature. This function 
	equals zero when tleaf does not change between iterations. """

	global lt

	tleaf_old = tleaf_in

	empirical_stomata(tleaf_old, mod, met_ts)
	tleaf_new = lt

	TleafFunc = tleaf_new - tleaf_old

	return TleafFunc

########################################################
#-------------------- Main code -----------------------#
########################################################

def main(in_file,met_data,start,end,model,g1):

	global met_ts, et, mod, g1_param

	mod = model
	g1_param = g1

	# initialize output data
	out_data = []
	print len(met_data)

	# loop through timesteps
	for i in range(len(met_data)):

		# get met data for timestep
		met_ts = met_timestep(met_data.iloc[i])

		# Run Stomatal Optimization Model 
		lt_out = optimize.brentq(TleafFunc, met_ts['temp']-20, met_ts['temp']+20)
		
		# prepare output data
		an = farquhar ( vcmax, vjmax, kc, ko, gamma1, resp, met_ts['par'], ci)
		agr = an+resp
		gsm = 1000*gs* met_ts['atmos_press'] / ( ( lt + 273.15 ) * Rcon ) # convert m/s to mmol/m2/s
		et = 18*et/1000 # convert from mol/m2/s to kg/m2/s
		
		# for comparison with SPA-WUEi scale by leaf area: la = 1/6 * LAI = .16667 * 3.4  
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
	#out_data.index = pd.date_range(start,end,freq='30min')

	return out_data

