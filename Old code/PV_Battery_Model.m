%   "An economic analysis of residential photovoltaic
%   systems with battery storage in the United States."

%   This model simulates residential PV/Battery systems over 30 years and
%   calculates electricity consumption, production, cost, etc.  A detailed
%   description can be found in the corresponding publication.  The model
%   was created with MATLAB R2016b.

%   Accompanying electricity price data, load data, and TMY3 data must be
%   placed in the same working directory as this script in order for the
%   model to run.  New data may be used as long as it is formatted to match
%   the original data files provided.

%   TMY3 files must be named 1.csv ... 50.csv for each state (currently alphabetical)
%   Corresponding load files must be named 1L.csv ... 51L.csv for each state

%   Last updated 2018-01-10

clearvars
clc
tic

nStates = 50;                           %number of states to run model for
nStart = 1;                             %state to begin with (numeric order based on file names)

%-----PV Parameters-----

etaPV_rated = 0.1534;                   %nameplate efficiency
costPV = 500;                           %avg cost of PV module w/o installation [$/kW]
capPV = 5;                              %installation capacity [kW]
areaPV = capPV / etaPV_rated;           %area of installed panels [m2] - based on 1 kW/m2 nominal irradiance
Gsc = 1367;                             %solar constant [W/m2]
gamma = 0;                              %surface azimuth angle [deg]
betaT = 0.00457;                        %coefficient for PV performance degradation with temp [1/C]
NOCT = 46.6;                            %parameter for calculating PV temp [C]
Tref = 25;                              %reference temperature for PV module [C]
etaMPP = 0.95;                          %maximum power point tracker efficiency
etaDust = 0.98;                         %efficiency reduction due to dust
etaDC = 0.98;                           %efficiency reduction due to DC losses
etaD = 0.008;                           %degradation ratio
etaI = 0.956;                           %inverter efficiency
costI = 160;                            %inverter cost [$/kW]
costIreplace = 120;                     %inverter replacement cost [$/kW]
replaceIyr = 15;                        %year to replace inverter
costITot = (costI * capPV);             %total inverter cost
costBOS = 500 * capPV;                  %balance of systems cost [$]
costInstall = (100+350+350+700)*capPV;  %installation cost [$]
costPermit = 100*capPV;                 %permitting cost [$]
costPVTot = costPV*capPV + costITot + costBOS + costInstall + costPermit;    %total cost of PV installation
costPVinstpW = costPVTot/capPV/1000;    %installed cost of PV [$/W]

%-----Storage Parameters-----

nomCapBat = 7;                          %nominal capacity of battery storage [kWh]
dod = 6.75/7;                           %allowed depth of discharge
etaStor = 0.90;                         %roundtrip storage efficiency
capStorRated = nomCapBat*dod;           %total initial storage capacity [kWh]
deltaCcal = 0.2*capStorRated/(15*365*24);    %hourly degradation due to calendric aging [kWh]
costBat = 392.86;                       %battery pack cost [$/kWh]
costBatreplace = costBat/2;             %cost to replace battery [$/kWh]
replaceBatYr = 15;                      %year in which to replace batteries
chargeMin = 0;                          %minimum storage charge (note that unallowed depth of...
                                            %discharge has been removed from total capacity) [kWh]
if nomCapBat > 0
    costStorInst = 1700;                %installation and BOS cost of batteries [$]
    costStorInstReplace = costStorInst*0.75;  %replacement installation and BOS cost of batteries [$]
else
    costStorInst = 0;
    costStorInstReplace = 0;
end
costStorTot = costBat*nomCapBat + costStorInst;   %total installed cost [$]


%-----Economic Parameters-----

nEc = 30;                               %period for economic analysis [yr]
m = 0.04;                               %loan interest rate
ITC = 0.3;                              %investment tax credit rate
dp = 0.2;                               %down payment ratio
tProp = 0.0;                            %property tax rate
mS = 0.01;                              %ratio insurance and maint. costs to investment
d = 0.06;                               %discount rate
tBar = 0.25;                            %effective income tax rate
r = 0.02;                               %inflation rate

invest = (costStorTot+costPVTot);       %total investment [$]
downPay = invest*dp;                    %down payment [$]

file = 'elecPrices.csv';                %electricity price projections file name
dataTable = readtable(file,'ReadVariableNames',0,'Delimiter',',','HeaderLines',3);  %read data [cents/kWh]
pElec = 0.01*table2array(dataTable(:,2:51));       %write electricity prices for each state to matrix [$/kWh]


%-----TMY3 Station Data & Load Data-----

hour = repmat((1:24).',365,1);                  %local time for each hour of year [hr]
n = reshape(repmat((1:365),24,1),8760,1);       %day of year for each hour of year [day]

TZ = NaN*ones(1,50);                            %initialize array variables
lat = NaN*ones(1,50);                           %latitude
long = NaN*ones(1,50);                          %longitude
beta = NaN*ones(1,50);
merid = NaN*ones(1,50);
I = NaN*ones(8760,50);                          %global irradiance
Id = NaN*ones(8760,50);                         %diffuse horizontal irradiance
T = NaN*ones(8760,50);                          %temperature
V = NaN*ones(8760,50);                          %velocity (wind speed)
rhoG = NaN*ones(8760,50);                       %albedo
eLoad = NaN*ones(8760,50);                      %hourly load

for i=1:nStates
    ii = i-1+nStart;                            %shifted index per starting state
    file = strcat(num2str(ii),'.csv');          %file name
    file = strcat('TMY3 Data/',file);           %adding folder name to point to right folder
    header = csvread(file,0,3,[0 3 0 5]);       %read header from csv file
    TZ(i) = header(1);                          %time zone
    lat(i) = header(2);                         %latitude [deg]
    long(i) = -header(3);                       %longitude in deg. west [deg]
    beta(i) = lat(i);                           %tilt angle from horizontal [deg]
    
    if (TZ(i)==-4), merid(i) = 60;              %set local meridian based on time zone [deg]
    elseif (TZ(i)==-5), merid(i) = 75;
    elseif (TZ(i)==-6), merid(i) = 90;
    elseif (TZ(i)==-7), merid(i) = 105;
    elseif (TZ(i)==-8), merid(i) = 120;
    elseif (TZ(i)==-9), merid(i) = 135;
    elseif (TZ(i)==-10), merid(i) = 150;
    end
    
    dataTable = readtable(file,'ReadVariableNames',0,'Delimiter',',','HeaderLines',2);
    I(:,i) = table2array(dataTable(:,5));       %global horizontal irradiance [W-hr/m2]
    Id(:,i) = table2array(dataTable(:,11));     %diffuse horizontal irradiance [W-hr/m2]
    T(:,i) = table2array(dataTable(:,32));      %temperature [C]
    V(:,i) = table2array(dataTable(:,47));      %wind speed [m/s]
    rhoG(:,i) = table2array(dataTable(:,62));   %albedo
    
    file = strcat(num2str(ii),'L.csv');         %load file name
    file = strcat('Load Data/',file);           %adding folder name to point to right folder
    dataTable = readtable(file,'ReadVariableNames',0,'Delimiter',',','HeaderLines',1);  %read file
    eLoad(:,i) = table2array(dataTable(:,2));   %hourly load [kWh]
end

disp('Completed processing input files.')
toc


%-----Run Simulation of PV/Battery System-----

eProdPV = NaN*ones(8760,50);                    %initialize energy produced by PV panels [kWh]
eProdPVTot = zeros*ones(nEc,50);                %initialize total energy produced by panels for year [kWh]
capStorHist = NaN*ones(nEc,50);                 %initialize yearly history of storage capacity [kWh]
eSysUt = zeros(nEc,50);                         %initialize total utilized energy produced by system [kWh]
eSysWa = zeros(nEc,50);                         %initialize total wasted energy produced by system [kWh]
eBatUse = zeros(nEc,50);                        %initialize total energy stored in batteries [kWh]
It = NaN*ones(8760,50);                         %initialize irradiance on tilted panel [W-hr/m2]
Tpv = NaN*ones(8760,60);                        %initialize temperature of PV module [C]
chargeHist = NaN*ones(8760,50);                 %initialize battery charge history [kWh]

for i=1:nStates

for j=1:8760

%--Solar Time and Hour Angles

B = (n(j)-1)*(360/365);                                 %factor for solar time [deg]
E = 229.2*(0.000075 + 0.001868*cosd(B) - 0.032077*sind(B) - 0.014615*cosd(2*B) - 0.04089*sind(2*B));   %factor for solar time
tSolDiff = 4*(merid(i)-long(i)) + E;                    %difference between standard and solar time [min]
tSol = hour(j) + tSolDiff/60;                           %solar time [hrs]
if tSol > 24
    tSol = (tSol-24);                                   %correct solar time if needed
elseif tSol < 0
    tSol = 24 + tSol;
end
omegaA = 15*(tSol - 1 - 12);                            %start of hour angle [deg]
omegaB = 15*(tSol - 12);                                %end of hour angle [deg]

delta = 23.45 * sind(360*(284+n(j))/365);               %solar declination angle
omegaSet = acosd(-tand(lat(i))*tand(delta));            %sunset angle [deg]
omegaRise = -omegaSet;                                  %sunrise angle [deg]

if omegaA < omegaRise && omegaB >= omegaRise
    omegaA = omegaRise;                                 %set start hour angle to sunrise if hour includes sunrise
end
if omegaA <= omegaSet && omegaB > omegaSet
    omegaB = omegaSet;                                  %set end hour angle to sunset if hour angle includes sunset
end

%--Rb Factor

aRb = (sind(delta)*sind(lat(i))*cosd(beta(i)) - sind(delta)*cosd(lat(i))*sind(beta(i))*cosd(gamma)) * (omegaB-omegaA)*pi/180 ...
	+ (cosd(delta)*cosd(lat(i))*cosd(beta(i)) + cosd(delta)*sind(lat(i))*sind(beta(i))*cosd(gamma)) * (sind(omegaB) - sind(omegaA)) ...
	- (cosd(delta)*sind(beta(i))*sind(gamma)) * (cosd(omegaB)-cosd(omegaA));

bRb = (cosd(lat(i))*cosd(delta))*(sind(omegaB)-sind(omegaA)) + (sind(lat(i))*sind(delta))*(omegaB-omegaA)*pi/180;
Rb = aRb/bRb;                                           %ratio of beam radiation on tilted surface to horizontal surface

if Rb<0                                                 %correct Rb when necessary
    Rb = 0;
elseif Rb>50
    Rb = 0;
end


%--Incident Irradiance on Tilted Panel

Io = (12*3600/pi)*Gsc*(1+0.033*cosd(360*n(j)/365)) ...
		* ((cosd(lat(i))*cosd(delta)*(sind(omegaB)-sind(omegaA))) ...
		+ pi*(omegaB-omegaA)*sind(lat(i))*sind(delta)/180);     %extraterrestrial radiation on horizontal surface [W-hr/m2]
    
Ib = I(j,i) - Id(j,i);                              %beam horizontal irradiance [W-hr/m2]
if Ib<0
    Ib=0;                                           %correct Ib for instances when TMY3 model overpredicts Id
end

Ai = Ib/Io;                                         %anisotropy index
f = sqrt(Ib/I(j,i));                                %modulating factor
if I(j,i)==0
    f=0;
end

It(j,i) = (Ib + Id(j,i)*Ai)*Rb ...
        + Id(j,i)*(1-Ai)*((1+cosd(beta(i)))/2)*(1+f*(sind(beta(i)/2))^3) ...
        + I(j,i)*rhoG(j,i)*((1-cosd(beta(i)))/2);   %irradiance on tilted plate [W-hr/m2]

if omegaB < omegaRise                               %set irradiance to zero at night
    It(j,i) = 0;
elseif omegaA > omegaSet
    It(j,i) = 0;
end
if Ai > 1 || Rb > 50
    It(j,1) = 0;
end
    

end
end

disp('Completed calculating incident irradiance on tilted panels.')
toc

capStor = capStorRated*ones(1,50);

for i=1:nStates
for k=1:nEc

if k==replaceBatYr
    capStor = capStorRated*ones(1,50);          %reset battery capacity when replaced
end
    
charge = chargeMin;                             %initialize battery charge [kWh]
chargeP = chargeMin;                            %initialize previous charge [kWh]
cycleMax = chargeMin;                           %initialize max charge on cycle [kWh]
cycleMin = chargeMin;                           %initialize min charge on cycle [kWh]
chargeDir = 0;                                  %initialize charge direction: -1 for discharge, 0 for steady, +1 for charge

for j=1:8760

%--Power Produced and Changes to Battery Storage

Tpv(j,i) = T(j,i) + It(j,i)*((NOCT-20)/800)...          %Temperature of PV module [C]
    *(1-etaPV_rated/0.9)*(9.5/(5.7+3.8*V(j,i)));
etaT = 1 - betaT*(Tpv(j,i) - Tref);                     %efficiency reduction due to PV temp
etaPV = etaPV_rated*etaT*etaDust*etaDC*etaMPP*(1-(k-1)*etaD);   %PV efficiency
eProdPV(j,i) = It(j,i)*etaPV*areaPV * (1/1000);         %electricity produced by PV panels in hour [kWh]
eProdPVTot(k,i) = eProdPVTot(k,i) + eProdPV(j,i);

if (eProdPV(j,i)*etaI <= eLoad(j,i))                    %if PV production is less than or equal to load
	eSysUt(k,i) = eSysUt(k,i) + eProdPV(j,i)*etaI;      %all PV production adds to total utilization
	loadEx = eLoad(j,i) - eProdPV(j,i)*etaI;            %excess load
    if (charge > chargeMin)                             %if battery is not at min charge
        if (loadEx >= (charge-chargeMin)*etaI*etaStor)  %and if excess load exceeds or is equal to available charge
	          eSysUt(k,i) = eSysUt(k,i) + (charge-chargeMin)*etaStor*etaI;      %avail. battery charge adds to total utilization
	          eBatUse(k,i) = eBatUse(k,i) + (charge-chargeMin)*etaStor*etaI;    %avail. battery chage adds to that used from battery
	          charge = chargeMin;                       %avail. battery charge goes to 0
         else                                           %else excess load is less than available charge
	          eSysUt(k,i) = eSysUt(k,i) + loadEx;       %excess load adds to total utilization
	          eBatUse(k,i) = eBatUse(k,i) + loadEx/(etaStor*etaI);     %excess load adds to that used from battery	
	          charge = charge - loadEx/(etaStor*etaI);  %charge decreases by amount equal to excess load
         end
    end
end

if (eProdPV(j,i)*etaI > eLoad(j,i))                     %if PV production is greater than load
	eSysUt(k,i) = eSysUt(k,i) + eLoad(j,i);             %system utilization increases by amount equal to load (used immediately)
	prodEx = eProdPV(j,i) - eLoad(j,i)/etaI;            %excess production
    if (prodEx > capStor(i)-charge)                     %if excess production exceeds uncharged capacity
	     eSysWa(k,i) = eSysWa(k,i) + (prodEx - (capStor(i)-charge));     %excess production above capacity is wasted
	     charge = capStor(i);                           %battery becomes fully charged
    else                                                %else excess production is less than or equal to uncharged capacity 
	     charge = charge + prodEx;                      %battery charge increases by excess production
    end
end
chargeHist(j,i) = charge;                               %store charge in history variable

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

capStor(i) = capStor(i) - deltaCcyc;                    %degrade capacity by cyclic aging
capStor(i) = capStor(i) - deltaCcal;                    %degrade capacity by calendric aging

if charge > capStor(i); charge = capStor(i); end

chargeP = charge;

end

capStorHist(k,i) = capStor(i);

end

end

SCR = eSysUt ./ eProdPVTot;                             %yearly self-consumption ratio
SSR = eSysUt ./ repmat(sum(eLoad),nEc,1);               %yearly self-sufficiency ratiio
BUR = eBatUse ./ (365*nomCapBat*ones(30,50));           %yearly battery utilization ratio
SCRtot = sum(eSysUt) ./ sum(eProdPVTot);                %total SCR
SSRtot = sum(eSysUt) ./ (nEc*sum(eLoad));               %total SSR
BURtot = sum(eBatUse) ./ sum(365*nomCapBat*ones(30,50));%total BUR


disp('Completed simulating PV/battery system for 30 years')
toc


%-----Financial Analysis of System-----
%%%%%%%%% Part of econAnalysis function %%%%%%%%%

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
end                                                             %**note that utilization/production is not truly discounted, this is an artifact of the LCOE eqn.

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

disp('Completed financial analysis')
toc


function p = PWF(N, i, d)
%   This function calculates the Present Worth Factor for N inflating costs
%   with inflation rate i and discount rate d.  The total present worth is
%   the annual cost multiplied by the Present Worth Factor.

if i==d
    p = N/(1+i);
else
    p = (1/(d-i)) * (1 - ((1+i)/(1+d))^N);
end

end