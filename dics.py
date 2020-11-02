#!/usr/bin/env python2.7

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
