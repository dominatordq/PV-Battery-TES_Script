function [LCOEsys, COEsys, LCCelecAvoid, LCOEgrid, LCOEcomp, LCOEpv, LCOEbat, Pbd] = econAnalysis(nEc,invest,downPay,ITC,m,tBar,d,mS,r,costIreplace,capPV,replaceIyr,costBatreplace,nomCapBat,costStorInstReplace,replaceBatYr,costPVTot,eSysUt,eProdPVTot,pElec)
%Calculation of economic figures of merit (LCOE of PV/battery/system, price for equivalent bi-directional metering, etc.)
%   Input: nEc = period for economic analysis [yr]
%   Input: invest = total investment [$]
%   Input: downPay = down payment [$]
%   Input: ITC = investment tax credit rate
%   Input: m = loan interest rate
%   Input: tBar = effective income tax rate
%   Input: d = discount rate
%   Input: mS = ratio insurance and maint. costs to investment
%   Input: r = inflation rate
%   Input: costIreplace = inverter replacement cost [$/kW]
%   Input: capPV = installation capacity [kW]
%   Input: replaceIyr = year to replace inverter
%   Input: costBatreplace = cost to replace battery [$/kWh]
%   Input: nomCapBat = nominal capacity of battery storage [kWh]
%   Input: costStorInstReplace = replacement installation and BOS cost of batteries [$]
%   Input: replaceBatYr = year in which to replace batteries
%   Input: costPVTot = total cost of PV installation
%   Input: eSysUt = total utilized energy produced by system [kWh]
%   Input: eProdPVTot = total energy produced by PV panels [kWh]
%   Input: pElec = electricity prices for each state [$/kWh]
%
%   Output: LCOEsys = LCOE from system
%   Output: COEsys = cost of electricity from system (non-discounted)
%   Output: LCCeleAvoid = life-cycle avoided cost of electricity [$]
%   Output: LCOEgrid = LCOE from grid
%   Output: LCOEcomp = LCOE comparison
%   Output: LCOEpv = production LCOE
%   Output: LCOEbat = battery LCOE
%   Output: Pbd = bi-directional sell-back price for battery parity

%-----Financial Analysis of System-----

P = NaN*ones(nEc+1,1);                                  %initialize remaining principal at beginning of year
taxSav = NaN*ones(nEc,1);                               %initialize tax savings each year
PWtaxSavY = NaN*ones(nEc,1);                            %initialize present worth of tax savings from each year

M = (invest - downPay - invest*ITC) / PWF(nEc,0,m);     %annual mortgage payment
P(1) = (invest - downPay - invest*ITC);                 %initial principal
mPay = NaN*ones(nEc,1);                                 %annual interest payment
for i = 1:nEc
    mPay(i) = P(i)*m;                                   %interest payment at end of year i
    P(i+1) = P(i) - (M-mPay(i));                        %calculate new remaining principal after payment
    taxSav(i) = mPay(i) * tBar;                         %tax savings each year
    PWtaxSavY(i) = taxSav(i) / (1+d)^i;                 %present worth of tax savings from each year
end

PWinvest = downPay - invest*ITC + M*PWF(nEc,0,d);       %present worth of investment
PWmisc = mS*invest*PWF(nEc,r,d);                        %present worth of insurance and misc costs
PWtaxSav = sum(PWtaxSavY);                              %total present worth of tax savings

PWreplaceI = costIreplace*capPV / (1+d)^replaceIyr;     %present worth of inverter replacement
PWreplaceBat = (costBatreplace*nomCapBat + costStorInstReplace) / (1+d)^replaceBatYr;    %present worth of battery replacement

PWsys = PWinvest + PWmisc - PWtaxSav + PWreplaceI + PWreplaceBat;       %net present worth of system costs
PWpv = (costPVTot/invest)*(PWinvest + PWmisc - PWtaxSav) + PWreplaceI;  %net present worth of PV portion only

eSysUtYPW = NaN*ones(nEc,50);                           %initialize present worth of yearly electricity utilization
eProdPVTotYPW = NaN*ones(nEc,50);                       %initialize present worth of yearly electricity production
for i = 1:nEc
    eSysUtYPW(i,:) = eSysUt(i,:)./(1+d)^i;              %present worth of yearly electricity utilization
    eProdPVTotYPW(i,:) = eProdPVTot(i,:)./(1+d)^i;      %present worth of yearly electricity production
end  

eSysUtPW = sum(eSysUtYPW,1);                            %calculate LCOE denominators
eProdPVTotPW = sum(eProdPVTotYPW,1);

LCOEsys = PWsys./eSysUtPW;                              %LCOE from system
COEsys = invest ./ sum(eSysUt);                         %cost of electricity from system (non-discounted)

PWelecAvoid = NaN*ones(nEc,50);                         %initialize present worth of avoided cost of electricity per year
for i = 1:nEc
    PWelecAvoid(i,:) = (eSysUt(i,:).*pElec(i,:))/(1+d)^i;    %present worth of avoided cost of electricity per year [$]
end
LCCelecAvoid = sum(PWelecAvoid);                        %life-cycle avoided cost of electricity [$]

LCOEgrid = LCCelecAvoid ./  eSysUtPW;                   %LCOE from grid
LCOEcomp = LCOEgrid - LCOEsys;                          %LCOE comparison

LCOEpv = PWpv./eProdPVTotPW;                            %production LCOE
LCOEbat = LCOEsys - LCOEpv;                             %battery LCOE
Pbd = LCOEgrid - LCOEbat;                               %bi-directional sell-back price for battery parity

end

