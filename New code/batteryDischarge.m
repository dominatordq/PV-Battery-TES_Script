function [chargeDirNew, currentCharge, cycleMin, cycleMax, deltaCcyc, capStorOut] = batteryDischarge(chargeDir,charge,chargeP,chargeMin,cycleMin,cycleMax,capStor,capStorRated,deltaCcal)
%This function will calculate updated discharge states of battery systems based on excess PV electricity, 
%charge states, technical parameters, and weather parameters.
%   Input: chargeDir = charge direction
%   Input: charge = battery charge [kWh]
%   Input: chargeP = previous charge [kWh]
%   Input: chargeMin = min charge on cycle [kWh]
%   Input: cycleMin = min charge on cycle [kWh]
%   Input: cycleMax = max charge on cycle [kWh]
%   Input: capStor = storage capacity [kWh]
%   Input: capStorRated = total initial storage capacity [kWh]
%   Input: deltaCcal = hourly degradation due to calendric aging [kWh]

if chargeDir == -1                                      %if battery was discharging
    if charge-chargeP > 0 || charge-chargeP == 0        %and if battery is now steady or charging
        cycleMin = chargeP;                             %set min cycle charge
        DOC = (cycleMax-cycleMin)/(capStorRated-chargeMin);     %calculate depth of cycle
        kCyc = real(4917.6*DOC^(-0.481));               %equivalent cycles until battery reaches 80% capacity
        deltaCcyc = 0.1*capStorRated/(kCyc);            %loss of capacity due to this half-cycle
        if charge-chargeP>0                             %update charge direction
            chargeDir = 1;
        else
            chargeDir = 0;
        end
    else
        deltaCcyc = 0;                                  %battery still discharging so no end to half-cycle
    end
elseif chargeDir == 1                                   %if battery was charging
    if charge-chargeP < 0 || charge-chargeP ==0         %and if battery is now steady or discharging
        cycleMax = chargeP;                             %set max cycle charge
        DOC = (cycleMax-cycleMin)/(capStorRated-chargeMin);     %calculate depth of cycle
        kCyc = real(4917.6*DOC^(-0.481));               %equivalent cycles until battery reaches 80% capacity
        deltaCcyc = 0.1*capStorRated/(kCyc);            %loss of capacity due to this half-cycle
        if charge-chargeP<0                             %update charge direction
            chargeDir = -1;
        else
            chargeDir = 0;
        end
    else
        deltaCcyc = 0;                                  %battery still charging so no end to half-cycle
    end
else                                                    %battery was steady
    deltaCcyc = 0;                                      %no end to half-cycle
    if charge-chargeP > 0                               %now charging
        cycleMin = chargeP;                             %set min cycle charge
        chargeDir = 1;                                  %update charge direction
    elseif charge-chargeP < 0                           %now discharging
        cycleMax = chargeP;                             %set max cycle charge
        chargeDir = -1;                                 %update charge direction
    end
end    

chargeDirNew = chargeDir;           %ouput the new charge direction

capStorOut = capStor - deltaCcyc;                    %degrade capacity by cyclic aging
capStorOut = capStorOut - deltaCcal;                 %degrade capacity by calendric aging
currentCharge = charge;
if charge > capStorOut; currentCharge = capStorOut; end

end

