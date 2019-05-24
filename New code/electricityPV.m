function [prodPV, prodPVTot, tempPV] = electricityPV(i,j,k,It,T,V,NOCT,areaPV,etaPV_rated,etaDust,etaDC,etaMPP,etaD,eProdPVTot,betaT,Tref)
%This function will calculate electricity produced by PV panels for each timestep, accounting for temperature, irradiance, and degradation effects. 
%Output will be total electricity produced for timestep
% Input: i, j, k = ith, jth, and kth index in for loop
% Input: It = irradiance on a tilted panel
% Input: T = temperature [C]
% Input: V = wind speed [m/s]
% Input: NOCT = parameter for calculating PV temp [C]
% Input: areaPV = area of installed panels [m2] - based on 1 kW/m2 nominal irradiance
% Input: etaPV_rated = nameplate efficiency
% Input: etaDust = efficiency reduction due to dust
% Input: etaDC = efficiency reduction due to DC losses
% Input: etaMPP = maximum power point tracker efficiency
% Input: etaD = degradation ratio
% Input: eProdPVTot = previous total energy produced by panels for year [kWh]
% Input: betaT = coefficient for PV performance degradation with temp [1/C]
% Input: Tref = reference temperature for PV module [C]
%
% Output: prodPV = electricity produced by PV panels in hour [kWh]
% Output: prodPVTot = total energy produced by panels for year [kWh]

tempPV = T(j,i) + It(j,i)*((NOCT-20)/800)...          %Temperature of PV module [C]
    *(1-etaPV_rated/0.9)*(9.5/(5.7+3.8*V(j,i)));
etaT = 1 - betaT*(tempPV - Tref);                    %efficiency reduction due to PV temp
etaPV = etaPV_rated*etaT*etaDust*etaDC*etaMPP*(1-(k-1)*etaD);   %PV efficiency
prodPV = It(j,i)*etaPV*areaPV * (1/1000);         %electricity produced by PV panels in hour [kWh]
prodPVTot = eProdPVTot(k,i) + prodPV;

end

