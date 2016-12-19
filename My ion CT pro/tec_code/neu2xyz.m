function XYZ = neu2xyz(NEU,lat,lon);

% NEU2XYZ	Converts local coordinates north, east, up (m) at location
%		lat, lon to ECEF coordinates X, Y, Z (m)
%		lat, lon in degrees
%		NEU =  input N,E,U matrix (n x 3), local coordinate system
%		XYZ = output X,Y,Z matrix (n x 3), 
%
%		XYZ = neu2xyz(NEU,lat,lon);
%

% convert lat, lon to radians
lat = lat .* (pi/180);
lon = lon .* (pi/180);

% build rotation matrix
clat = cos(lat); slat = sin(lat);
clon = cos(lon); slon = sin(lon);
R = [ -slat.*clon  -slon  clat.*clon;
      -slat.*slon   clon  clat.*slon;
       clat         0     slat];

% apply rotation
XYZ = R * NEU';

