function [SCR, SSR, BUR, SCRtot, SSRtot, BURtot] = techAnalysis(eSysUt,eBatUse,eProdPVTot,eLoad,nomCapBat,nEc)
%Calculation of technical figures of merit (self-consumption ratio, battery utilization ratio, etc.) 
%   Input: eSysUt = utilized energy produced by system [kWh]
%   Input: eBatUse = energy stored in batteries [kWh]
%   Input: eProdPVTot = total energy produced by PV panels [kWh]
%   Input: eLoad = hourly load
%   Input: nomCapBat = nominal capacity of battery storage [kWh]
%   Input: nEc = period for economic analysis [yr]
%
%   Output: SCR = yearly self-consumption ratio
%   Output: SSR = yearly self-sufficiency ratiio
%   Output: BUR = yearly battery utilization ratio
%   Output: SCRtot = total SCR
%   Output: SSRtot = total SSR
%   Output: BURtot = total BUR

SCR = eSysUt ./ eProdPVTot;                             %yearly self-consumption ratio
SSR = eSysUt ./ repmat(sum(eLoad),nEc,1);               %yearly self-sufficiency ratiio
BUR = eBatUse ./ (365*nomCapBat*ones(30,50));           %yearly battery utilization ratio
SCRtot = sum(eSysUt) ./ sum(eProdPVTot);                %total SCR
SSRtot = sum(eSysUt) ./ (nEc*sum(eLoad));               %total SSR
BURtot = sum(eBatUse) ./ sum(365*nomCapBat*ones(30,50));%total BUR

end

