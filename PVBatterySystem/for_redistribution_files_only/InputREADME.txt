------------Legend--------------
W = watt
kW = kilowatt
kWh = kilowatt hour(s)
$ = US dollars
m^2 = meters squared
deg = degrees
C = Celsius 
yr = year
--------------------------------


nStates: number of states to run model for

nStart: state to begin with (numeric order based on file names)


-------------PV Parameters-------------

etaPV_rated: nameplate efficiency

costPV: average cost of PV module without installation [$/kW]

capPV: installation capacity [kW]

areaPV: area of installed panels [m^2] - based on 1 kW/m^2 nominal irradiance (capPV / etaPV_rated)

Gsc: solar constant [W/m^2]

gamma: surface azimuth angle [deg]

betaT: coefficient for PV performance degradation with temp [1/C]

NOCT: parameter for calculating PV temp [C]

Tref: reference temperature for PV module [C]

etaMPP: maximum power point tracker efficiency

etaDust: efficiency reduction due to dust

etaDC: efficiency reduction due to DC losses

etaD: degradation ratio

etaI: inverter efficiency

costI: inverter cost [$/kW]

costIreplace: inverter replacement cost [$/kW]

replaceIyr: year to replace inverter

costITot: total inverter cost

costBOS: balance of systems cost (500 * capPV) [$]

costInstall: installation cost ((100+350+350+700) * capPV) [$]

costPermit: permitting cost (100 * capPV) [$]

costPVTot: total cost of PV installation (costPV*capPV + costITot + costBOS + costInstall + costPermit)

costPVinstpW: installed cost of PV (costPVTot/capPV/1000) [$]


-------------Storage Parameters-------------

nomCapBat: nominal capacity of battery storage [kWh]

dod: allowed depth of discharge

etaStor: roundtrip storage efficiency

capStorRated: total initial storage capacity (nomCapBat*dod) [kWh]

deltaCcal: hourly degradation due to calendric aging (0.2 * capStorRated / * (15*365*24)) [kWh]

costBat: battery pack cost [$/kWh]

costBatreplace: cost to replace battery (costBat / 2) [$/kWh]

replaceBatYr: year in which to replace batteries

chargeMin: minimum storage charge


-------------Economic Parameters-------------

nEc: period for economic analysis [yr] 

m: loan interest rate

ITC: investment tax credit rate

dp: down payment ratio

tProp: property tax rate

mS: ratio insurance and maint. costs to investment

d: discount rate

tBar: effective income tax rate

r: inflation rate
