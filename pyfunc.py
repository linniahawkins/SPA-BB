#!/usr/bin/env python2.7

import netCDF4 as nc
import numpy as np
import pandas as pd
import os

def farquhar( vcmax, vjmax, kc, ko, gamma1, resp, par, ci ):
	"""Calculate leaf internal co2 concentration (ci) from stomatal conductance (gs)
	assuming a leaf temperature of 25C

    Parameters : vcmax   : stomatal conductance (mol H2O/m2 leaf/s)
    			 vjmax   : ambient CO2 concentration (umol/mol)
    			 kc      : Michaelis-Menten constants for CO2
    			 ko      : and oxygen
    			 par     : PAR
    			 resp    : respiration
    			 gamma1  : CO2 compensation point
    			 ci      : internal co2 concentration (umol/mol)
                 
    Returns    : anx     : metabolic net photosynthetic rate (umolC/m2 ground area /s)"""

	oi = 210 # o2 partial pressure umol/mol
	alphaj = 0.385 # initial slope of quantum response curve
	theta = 0.7 # curvature of quantum response curve

	# determine Rubisco limited carboxylation rate
	wc = ( vcmax*ci ) / (ci + kc * ( 1+oi / ko ))

	# determine potential rate of RuBP regeneration
	bee = alphaj * par + vjmax
	vj = ( bee - np.sqrt( bee**2 - 4 * theta * alphaj * par * vjmax )) / ( 2*theta )

	# determine RuBP regeneration limited carboxylation rate
	wj = vj * ci / ( 4.5 * ci + 10.5*gamma1)

	# determine limiting rate, carboxylation or regeneration
	vc = min( wc, wj )

	# net photosynthetic rate, carboxylation regeneration
	farquhar = ( vc * ( 1 - gamma1 / ci ) - resp )

	return farquhar


def leaf_temperature( gs, netrad, temp, wdef, gbb):
    """ determines leaf temperature from stomatal conductance

    Parameters  :   gs   # stomatal conductance (m/s)
       				netrad # canopy level net radiation (kW.m-2)
    				temp # air temperature in C
    				wdef # water deficit of air (kg.m-3)
    				gbb # canopy level boundary conductance for water vapour (m.s-1) 
    Returns		:	leaf temperature (C) """

    # set constants
    cp_air = 1004 # Specific heat capacity of air; J.kg-1.K-1
    air_density_kg = 353 / (temp + 273.15 )
    gbh = .03  #! boundary layer conductance for heat for a given canopy layer (m.s-1)
    boltz_kW = 5.67*10**-11 
    emiss = 0.959

    ta = temp + 273.15 # convert to K
    rho =  air_density_kg * 10**3 # density of air g/m3
    cp = cp_air*10**-6 # unit conversion

    # slope of saturation vapor pressure curve (t-dependent)
    s = 6.1078 * 17.269 * 237.3 * np.exp( 17.269*temp / (237.3+temp))
    slope = 0.1 * ( s / (237.3+temp )**2)
    psych = 0.1 * ( 0.646 * np.exp( 0.00097*temp ))
    lambd = 0.001 * (2501 - 2.364*temp)

    # convert water deficit to vapor presure deficit
    de = wdef * lambd*psych / (rho*cp)
    gr = 4*emiss * boltz_kW * (ta**3) / (rho * cp )
    ghr = gr + 2 * gbh
    rhr = 1 / ghr
    rt = 1/gs + 1/gbb
    denom = psych * rt + slope*rhr
    diff = rhr*rt*psych*netrad / (rho*cp*denom) - rhr*de / denom

    leaf_temperature = diff+temp

    return leaf_temperature


def evap( gs, lt, netrad, wdef, gbb ):
    """determine evapotranspiration rate (g m-2 s-1) 
    from q (kW m-2), lt (oC), wdef (g m-3) and gbb and gs (m s-1)    

    Parameters  :    gbb  # canopy level boundary conductance for water vapour (m.s-1)
                     gs   # incoming stomatal conductance (m.s-1)
                     netrad    # canopy level net radiation (kW.m-2)
                     lt   # incoming leaf temperature (oC)
                     wdef # water deficit of air (kg.m-3)"""

    # slope of saturation vapour pressure curve (t-dependent)
    s      = 6.1078 * 17.269 * 237.3 * np.exp( 17.269 * lt / ( 237.3 + lt ) )
    slope  = 0.1 * ( s / ( 237.3 + lt )**2 )      # (kPa K-1)
    psych  = 0.0646 * np.exp( 0.00097 * lt )         # psych is temp-dependent (kPa K-1)
    eps    = slope / psych                            # response of epsilon to temp
    lambd = ( 2.5010 - 0.002364 * lt )           # latent heat of vapourisation (KJ g-1)
    evap   = ( eps * netrad / lambd + wdef * gbb ) / ( 1 + eps + gbb / gs ) # (g m-2 s-1)

    return evap


def diffusion( gs, ci, et , gbb , lt, atmos_press, co2 ):
	"""diffusion limited assimilation rate (umol.C.m-2 ground area s-1)
	Parameters  :   gs   # incoming stomatal conductance (m.s-1)
					ci   # internal CO2 concentration (umol/mol)
					et   # evapotranspiration 
					gbb  # boundary layer conductance (m/s)
					lt   # leaf temperature (C)
					atmos_press # atmospheric pressure(Pa)
					co2  # ambient CO2 (umol/mol)
	Returns		:	adx diffusion based photosynthetic rate (umolC/m2 ground area /s) """

	gi = 1 # mesophyll conductance (m/s)
	Rcon = 8.3144 # universal gas constant (J/mol/K)

	# total leaf conductance (converts from m/s to mol/s)
	convert = atmos_press / ( Rcon * (273.15 + lt))
	gt = convert / (1.65/gs + 1.37/gbb + 1/gi)

	diffusion = ( gt - 0.5*et ) * co2 - ( gt + 0.5*et ) * ci

	return diffusion


def quadratic(a,b,c):
	
	if (b < 0 ):
		q = -0.5 * (b - np.sqrt(b*b - 4*a*c))
	else:
	    q = -0.5 * (b + np.sqrt(b*b - 4*a*c))

	r1 = q/a

	if (q == 0):
		r2 = 1*(10**36)
	else:
		r2 = c / q

	return [r1,r2]


def arrhenious( a, b, lt):
	""" Temperature adjustment for Michaelis-Menten coeefficients for CO2 (kc)
	and 02 (k0) and co2 compensation point
			a  	: saturation
			b 	: 	half saturation concentration
			t 	: 	leaf temperature """

	arrhenious = a * np.exp( b * (lt - 25) / (lt + 273.15))

	return arrhenious


def tempmet( max_temp , opt , q , lt):
	"""Apply non-gaussian temperature modifications on carboxylation and 
	electron-transport rates
		max_temp 	:	max temperature
		opt 	 	: 	metabolic optimum temperature
		q 			: 	kurtosis
		lt 			: 	leaf temperature (C)"""

	if ( lt > max_temp ): 
		tempmet = 0
	else:
		dummy = (max_temp - lt ) / (max_temp - opt )
		dummy = np.exp( np.log( dummy ) * q * ( max_temp - opt ))
		tempmet = dummy * np.exp(q * (lt-opt))

	return tempmet
	


def boundary( temp, tower_ht, atmos_press, wind_spd, dimen ):
	""" Calculate boundary layer conductances: layer height hard coded to 14m 
	Parameters	temp 		:	atmospheric temperature (C)
				tower_ht 	: 	tower height (m)
				atmos_press : 	atmospheric pressure (Pa)
				wind_spd	: wind speed (m/s)
				dimen		: leaf characteristic dimension 
	Returns 	boundary layer conductance to water vapor (m/s) and heat"""

	layer_ht = 14
	alpha  = 4

	Dwv   = 0.0000242 * ( ( ( temp + 273.15 ) / 293.15 )**1.75 ) * 101300 / atmos_press
	mult   = np.exp( alpha * ( layer_ht / tower_ht - 1 ) )
	u = wind_spd * mult
	thick  = 0.004 * ( dimen / u )**0.5
	# conductance to water vapour (m s-1 - not mol m-2 s-1) i.e. remove P/RT
	gbw = Dwv / thick
	# approximate boundary layer conductance to heat
	leaf_heat_conductance = gbw * 0.93

	return gbw, leaf_heat_conductance











