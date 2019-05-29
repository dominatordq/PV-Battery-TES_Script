function [irradiance, isTrue] = irradiancePV(i,j,I,Id,rhoG,beta,gamma,merid,lat,long,Gsc,n)
%This function will process solar/weather data to determine irradiance on tilted PV panel.
% Input: i = ith index (row)
% Input: j = jth index (col)
% Input: I = global irradiance
% Input: Id = diffuse horizontal irradiance
% Input: rhoG = albedo
% Input: beta = tilt angle
% Input: gamma
% Input: merid = meridian
% Input: lat = latitude
% Input: long = longitude
% Input: Gsc = solar constant [W/m2]
% Input: n = day of year for each hour of year [day]
%
% Output: irradiance = irradiance on a tilted panel
% Output: isTrue = 1 or 0 depending if Ai > 1 or Rb > 50

hour = repmat((1:24).',365,1);                  %local time for each hour of year [hr]
isTrue = 0;                                     %1 if Ai > 1 || Rb > 50

%--Solar Time and Hour Angles

B = (n-1)*(360/365);                                 %factor for solar time [deg]
E = 229.2*(0.000075 + 0.001868*cosd(B) - 0.032077*sind(B) - 0.014615*cosd(2*B) - 0.04089*sind(2*B));   %factor for solar time
tSolDiff = 4*(merid-long) + E;                    %difference between standard and solar time [min]
tSol = hour(j) + tSolDiff/60;                           %solar time [hrs]
if tSol > 24
    tSol = (tSol-24);                                   %correct solar time if needed
elseif tSol < 0
    tSol = 24 + tSol;
end
omegaA = 15*(tSol - 1 - 12);                            %start of hour angle [deg]
omegaB = 15*(tSol - 12);                                %end of hour angle [deg]

delta = 23.45 * sind(360*(284+n)/365);               %solar declination angle
omegaSet = acosd(-tand(lat)*tand(delta));            %sunset angle [deg]
omegaRise = -omegaSet;                                  %sunrise angle [deg]

if omegaA < omegaRise && omegaB >= omegaRise
    omegaA = omegaRise;                                 %set start hour angle to sunrise if hour includes sunrise
end
if omegaA <= omegaSet && omegaB > omegaSet
    omegaB = omegaSet;                                  %set end hour angle to sunset if hour angle includes sunset
end

%--Rb Factor

aRb = (sind(delta)*sind(lat)*cosd(beta) - sind(delta)*cosd(lat)*sind(beta)*cosd(gamma)) * (omegaB-omegaA)*pi/180 ...
	+ (cosd(delta)*cosd(lat)*cosd(beta) + cosd(delta)*sind(lat)*sind(beta)*cosd(gamma)) * (sind(omegaB) - sind(omegaA)) ...
	- (cosd(delta)*sind(beta)*sind(gamma)) * (cosd(omegaB)-cosd(omegaA));

bRb = (cosd(lat)*cosd(delta))*(sind(omegaB)-sind(omegaA)) + (sind(lat)*sind(delta))*(omegaB-omegaA)*pi/180;
Rb = aRb/bRb;                                           %ratio of beam radiation on tilted surface to horizontal surface

if Rb<0                                                 %correct Rb when necessary
    Rb = 0;
elseif Rb>50
    Rb = 0;
end

%--Incident Irradiance on Tilted Panel

Io = (12*3600/pi)*Gsc*(1+0.033*cosd(360*n/365)) ...
		* ((cosd(lat)*cosd(delta)*(sind(omegaB)-sind(omegaA))) ...
		+ pi*(omegaB-omegaA)*sind(lat)*sind(delta)/180);     %extraterrestrial radiation on horizontal surface [W-hr/m2]
    
Ib = I - Id;                              %beam horizontal irradiance [W-hr/m2]
if Ib<0
    Ib=0;                                           %correct Ib for instances when TMY3 model overpredicts Id
end

Ai = Ib/Io;                                         %anisotropy index
f = sqrt(Ib/I);                                %modulating factor
if I==0
    f=0;
end

irradiance = (Ib + Id*Ai)*Rb ...
        + Id*(1-Ai)*((1+cosd(beta))/2)*(1+f*(sind(beta/2))^3) ...
        + I*rhoG*((1-cosd(beta))/2);   %irradiance on tilted plate [W-hr/m2]

if omegaB < omegaRise                               %set irradiance to zero at night
    irradiance = 0;
elseif omegaA > omegaSet
    irradiance = 0;
end
if Ai > 1 || Rb > 50
    irradiance = 0;
    isTrue = 1;
end

end

