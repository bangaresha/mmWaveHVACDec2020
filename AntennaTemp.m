%%
% PhiMax = 2*pi; % Max Azimuth Angle
% ThetaMax = -pi/2;  % Max Elevation Angle

function [antennaGainAE] = AntennaTemp (antennaGainRes,demoMode)
Length = 1;
% phi= linspace(0,pi/2,antennaGainRes); %pi = -28.5
% theta= linspace(0,pi,antennaGainRes);
phi= linspace(0,2*pi,antennaGainRes);
theta= linspace(0,pi,antennaGainRes);
[THETA,PHI]=meshgrid(theta,phi);
R=abs((cos(pi.*Length.*cos(THETA))- cos(pi.*Length))./(sin(THETA+eps)));
x=R.*sin(THETA).*cos(PHI);
y=R.*sin(THETA).*sin(PHI);
z=R.*cos(THETA);

if demoMode == 1

    figure
    mesh(x,y,z);
end
antennaGainAE = R;
title(['Antenna Pattern, Resolution = ',num2str(antennaGainRes)]);