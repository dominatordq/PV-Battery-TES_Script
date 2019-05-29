function [sysUt, sysWa, batUse, currentCharge] = batteryCharge(eProdPV,eLoad,charge,chargeMin,etaI,etaStor,eSysU,eSysW,eBatU,capStor)
%This function will calculate updated charge states of battery systems based on excess PV electricity, 
%charge states, technical parameters, and weather parameters.
%   Input: i, j = ith, jth, indices
%   Input: eProdPV = energy produced by PV panels [kWh]
%   Input: eLoad = hourly load
%   Input: charge = battery charge [kWh]
%   Input: chargeMin = min charge on cycle [kWh]
%   Input: etaI = inverter efficiency
%   Input: etaStor = roundtrip storage efficiency
%   Input: eSysUt = total utilized energy produced by system [kWh]
%   Input: eSysWa = total wasted energy produced by system [kWh]
%   Input: eBatUse = total energy stored in batteries [kWh]
%   Input: capStor = capacity storage
%   Input: nStates = number of states
%   Input: nEc = period for economic analysis [yr]
%   
%   Output: sysUt = utilized energy produced by system [kWh]
%   Output: sysWa = wasted energy produced by system [kWh]
%   Output: batUse = energy stored in batteries [kWh]
%   Output: currentCharge = current battery charge [kWh]

% nStates = 50; 
% nEc = 30;
% eSysUt(k,i) = 0;
% eBatUse(k,i) = 0;
eSysW = 0;      %initialize wasted energy just in case it never gets assigned
eBatU = 0;      %initialize battery use just in case it never gets assigned

if (eProdPV*etaI <= eLoad)                    %if PV production is less than or equal to load
    eSysU = eProdPV*etaI;      %all PV production adds to total utilization
    loadEx = eLoad - eProdPV*etaI;            %excess load
    if (charge > chargeMin)                             %if battery is not at min charge
        if (loadEx >= (charge-chargeMin)*etaI*etaStor)  %and if excess load exceeds or is equal to available charge
              eSysU = eSysU +(charge-chargeMin)*etaStor*etaI;      %avail. battery charge adds to total utilization
	          eBatU = (charge-chargeMin)*etaStor*etaI;  %avail. battery charge adds to that used from battery
              charge = chargeMin;                       %avail. battery charge goes to 0
         else                                           %else excess load is less than available charge
              eSysU = eSysU + loadEx;       %excess load adds to total utilization
	          eBatU = loadEx/(etaStor*etaI);     %excess load adds to that used from battery	
	          charge = charge - loadEx/(etaStor*etaI);  %charge decreases by amount equal to excess load
         end
    end
end

if (eProdPV*etaI > eLoad)                     %if PV production is greater than load
	eSysU = eLoad;             %system utilization increases by amount equal to load (used immediately)
    prodEx = eProdPV - eLoad/etaI;            %excess production
    if (prodEx > capStor-charge)                     %if excess production exceeds uncharged capacity
	     eSysW = (prodEx - (capStor-charge));        %excess production above capacity is wasted
	     charge = capStor;                           %battery becomes fully charged
    else                                                %else excess production is less than or equal to uncharged capacity 
	     charge = charge + prodEx;                      %battery charge increases by excess production
    end
end

sysUt = eSysU;
sysWa = eSysW;
batUse = eBatU;
currentCharge = charge;


end
