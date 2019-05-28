function [sysUt,batUse] = batteryCharge(i,j,k,charge,chargeMin,loadEx,etaI,etaStor,eSysUt,eBatUse)
%This function will calculate updated charge states of battery systems based on excess PV electricity, 
%charge states, technical parameters, and weather parameters.
%   Input: i, j, k = ith, jth, and kth indices
%   Input: charge = battery charge [kWh]
%   Input: chargeMin = min charge on cycle [kWh]
%   Input: loadEx = excess load
%   Input: etaI = inverter efficiency
%   Input: etaStor = roundtrip storage efficiency
%   
%   Output: sysUt = utilized energy produced by system [kWh]
%   Output: batUse = total energy stored in batteries [kWh]



end

