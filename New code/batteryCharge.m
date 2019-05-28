function [sysUt, sysWa, batUse, currentCharge] = batteryCharge(i,j,k,eProdPV,eLoad,charge,chargeMin,chargeHist,etaI,etaStor,eSysUt,eSysWa,eBatUse,capStor)
%This function will calculate updated charge states of battery systems based on excess PV electricity, 
%charge states, technical parameters, and weather parameters.
%   Input: i, j, k = ith, jth, and kth indices
%   Input: eProdPV = energy produced by PV panels [kWh]
%   Input: eLoad = hourly load
%   Input: charge = battery charge [kWh]
%   Input: chargeMin = min charge on cycle [kWh]
%   Input: chargeHist = battery charge history [kWh]
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
%   Output: chargeHistOut = updated battery charge history [kWh]

% nStates = 50; 
% nEc = 30;
batUse = 0;
sysWa = 0;
if (eProdPV(j,i)*etaI <= eLoad(j,i))                    %if PV production is less than or equal to load
	sysUt = eSysUt(k,i) + eProdPV(j,i)*etaI;      %all PV production adds to total utilization
	loadEx = eLoad(j,i) - eProdPV(j,i)*etaI;            %excess load
    if (charge > chargeMin)                             %if battery is not at min charge
        if (loadEx >= (charge-chargeMin)*etaI*etaStor)  %and if excess load exceeds or is equal to available charge
	          sysUt = sysUt + eSysUt(k,i) + (charge-chargeMin)*etaStor*etaI;      %avail. battery charge adds to total utilization
	          batUse = eBatUse(k,i) + (charge-chargeMin)*etaStor*etaI;    %avail. battery charge adds to that used from battery
	          charge = chargeMin;                       %avail. battery charge goes to 0
         else                                           %else excess load is less than available charge
	          sysUt = eSysUt(k,i) + loadEx;       %excess load adds to total utilization
	          batUse = eBatUse(k,i) + loadEx/(etaStor*etaI);     %excess load adds to that used from battery	
	          charge = charge - loadEx/(etaStor*etaI);  %charge decreases by amount equal to excess load
         end
    end
end

if (eProdPV(j,i)*etaI > eLoad(j,i))                     %if PV production is greater than load
	sysUt = eSysUt(k,i) + eLoad(j,i);             %system utilization increases by amount equal to load (used immediately)
	prodEx = eProdPV(j,i) - eLoad(j,i)/etaI;            %excess production
    if (prodEx > capStor(i)-charge)                     %if excess production exceeds uncharged capacity
	     sysWa(k,i) = eSysWa(k,i) + (prodEx - (capStor(i)-charge));     %excess production above capacity is wasted
	     charge = capStor(i);                           %battery becomes fully charged
    else                                                %else excess production is less than or equal to uncharged capacity 
	     charge = charge + prodEx;                      %battery charge increases by excess production
    end
end
currentCharge = charge;

end
