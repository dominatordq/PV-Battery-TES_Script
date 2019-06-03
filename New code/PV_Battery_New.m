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
%   Corresponding load files must be named 1L.csv ... 50L.csv for each state

%   Last updated 2019-06-03

clearvars
clc

filename = '';
if (isempty(filename))
    prompt = 'Enter the name of your input file (ex: Inputs.xlsx): ';   %this if statement prevents program from looping
end

filename = input(prompt,'s');               %save input into a variable
inputTable = xlsread(num2str(filename));    %read in excel file with user inputs

tic

nStates = inputTable(1);                %number of states to run model for
nStart = inputTable(2);                 %state to begin with (numeric order based on file names)

%-----PV Parameters-----

etaPV_rated = inputTable(3);            %nameplate efficiency
costPV = inputTable(4);                 %avg cost of PV module w/o installation [$/kW]
capPV = inputTable(5);                  %installation capacity [kW]
areaPV = inputTable(6);                 %area of installed panels [m2] - based on 1 kW/m2 nominal irradiance
Gsc = inputTable(7);                    %solar constant [W/m2]
gamma = inputTable(8);                  %surface azimuth angle [deg]
betaT = inputTable(9);                  %coefficient for PV performance degradation with temp [1/C]
NOCT = inputTable(10);                  %parameter for calculating PV temp [C]
Tref = inputTable(11);                  %reference temperature for PV module [C]
etaMPP = inputTable(12);                %maximum power point tracker efficiency
etaDust = inputTable(13);               %efficiency reduction due to dust
etaDC = inputTable(14);                 %efficiency reduction due to DC losses
etaD = inputTable(15);                  %degradation ratio
etaI = inputTable(16);                  %inverter efficiency
costI = inputTable(17);                 %inverter cost [$/kW]
costIreplace = inputTable(18);          %inverter replacement cost [$/kW]
replaceIyr = inputTable(19);            %year to replace inverter
costITot = inputTable(20);              %total inverter cost
costBOS = inputTable(21);               %balance of systems cost [$]
costInstall = inputTable(22);           %installation cost [$]
costPermit = inputTable(23);            %permitting cost [$]
costPVTot = inputTable(24);             %total cost of PV installation
costPVinstpW = inputTable(25);          %installed cost of PV [$/W]

%-----Storage Parameters-----

nomCapBat = inputTable(26);             %nominal capacity of battery storage [kWh]
dod = inputTable(27);                   %allowed depth of discharge
etaStor = inputTable(28);               %roundtrip storage efficiency
capStorRated = inputTable(29);          %total initial storage capacity [kWh]
deltaCcal = inputTable(30);             %hourly degradation due to calendric aging [kWh]
costBat = inputTable(31);               %battery pack cost [$/kWh]
costBatreplace = inputTable(32);        %cost to replace battery [$/kWh]
replaceBatYr = inputTable(33);          %year in which to replace batteries
chargeMin = inputTable(34);             %minimum storage charge (note that unallowed depth of...
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

nEc = inputTable(35);                   %period for economic analysis [yr]
m = inputTable(36);                     %loan interest rate
ITC = inputTable(37);                   %investment tax credit rate
dp = inputTable(38);                    %down payment ratio
tProp = inputTable(39);                 %property tax rate
mS = inputTable(40);                    %ratio insurance and maint. costs to investment
d = inputTable(41);                     %discount rate
tBar = inputTable(42);                  %effective income tax rate
r = inputTable(43);                     %inflation rate
                                            
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
    file = strcat(num2str(ii),'.CSV');          %file name
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

[It] = irradiancePV(I,Id,It,rhoG,beta,gamma,merid,lat,long,Gsc,n,nStates);  %call irradiancePV

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
%call electrivityPV to update efficiency/electricity 
[prodPV, prodPVTot, tempPV] = electricityPV(k,It(j,i),T(j,i),V(j,i),NOCT,areaPV,etaPV_rated,etaDust,etaDC,etaMPP,etaD,eProdPVTot(k,i),betaT,Tref); %call electricityPV
%set the outputs to their corresponding vectors
eProdPV(j,i) = prodPV;
eProdPVTot(k,i) = prodPVTot;
Tpv(j,i) = tempPV;
                 
[sysUt, sysWa, batUse, currentCharge] = batteryCharge(eProdPV(j,i),eLoad(j,i),charge,chargeMin,etaI,etaStor,capStor(i)); %call batteryCharge  
eSysUt(k,i) = eSysUt(k,i) + sysUt;      %add the sysUt output to the eSysUt array/matrix if there is a change
eSysWa(k,i) = eSysWa(k,i) + sysWa;      %add the sysWa output to the eSysWa array/matrix if there is a change
eBatUse(k,i) = eBatUse(k,i) + batUse;   %add the batUse output to the eBatUse array/matrix if there is a change
charge = currentCharge;         %store charge output 

chargeHist(j,i) = charge;    %store charge in history variable

[chargeDirNew, currentCharge, cycleMin, cycleMax, capStorOut] = batteryDischarge(chargeDir,charge,chargeP,chargeMin,cycleMin,cycleMax,capStor(i),capStorRated,deltaCcal); %call batteryDischarge
chargeDir = chargeDirNew;   %update charge direction
capStor(i) = capStorOut;    %update capacity storage array
charge = currentCharge;     %update charge
chargeP = charge;           %update previous charge

end
capStorHist(k,i) = capStor(i);  %update storage capacity history
end
end

[SCR, SSR, BUR, SCRtot, SSRtot, BURtot] = techAnalysis(eSysUt,eBatUse,eProdPVTot,eLoad,nomCapBat,nEc);  %call techAnalysis

disp('Completed simulating PV/battery system for 30 years')
toc

[LCOEsys, COEsys, LCCelecAvoid, LCOEgrid, LCOEcomp, LCOEpv, LCOEbat, Pbd] = econAnalysis(nEc,invest,downPay,ITC,m,tBar,d,mS,r,costIreplace,capPV,replaceIyr,costBatreplace,nomCapBat,costStorInstReplace,replaceBatYr,costPVTot,eSysUt,eProdPVTot,pElec); %call econAnalysis

disp('Completed financial analysis')
toc