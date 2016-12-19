function XYZ = elleza(AZ,EL,LE,P);

% ELLEZA	Compute XYZ ECEF coordinates (in m) of vector defined
%		by local azimuth (AZ), elevation (EL), and length (LE)
%		at point P (= [X Y Z], ECEF, m)
%		AZ, EL in degres
%		LE in meters
%
%		XYZ = elleza(AZ,EL,LE,P);
%

% convert elevation and azimuth to radians
EL = EL .* (pi/180);
AZ = AZ .* (pi/180);

% convert elevation, azimuth, length to local NEU coordinates in m
n = cos(EL) .* cos(AZ) .* LE;
e = cos(EL) .* sin(AZ) .* LE;
u = sin(EL) .* LE;

% get latitude and longitude of reference point
E = xyz2wgs([0 P(1) P(2) P(3)]);
lon = E(2);
lat = E(3);

% convert local NEU coordinates to ECEF
XYZ = neu2xyz([n e u],lat,lon);
