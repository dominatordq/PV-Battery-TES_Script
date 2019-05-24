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

%   Last updated 2019-05-23

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
It = NaN*ones(8760,50);                         %irradiance on a tilted panel
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

[irr, isTrue] = irradiancePV(i,j,I,Id,rhoG,beta,gamma,merid,lat,long,Gsc,n);  %call irradiancePV function
It(j,i) = irr;
if (isTrue == 1)    %if Ai > 1 || Rb > 50
    It(j,1) = 0;
end

end
end

disp('Completed calculating incident irradiance on tilted panels.')
toc

