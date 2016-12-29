time = datenum([2007 7 17 6 30 0]);
lat_start = 42.17; % Geodetic latitudes in degrees.
lon_start = 128.0; % Geodetic longitudes in degrees.
Re=6371.2; %km
ALTITUDE =205:5:400; % Altitude in km.

load igrfcoefs.mat
[BX, BY, BZ]=igrf(time, lat_start, lon_start, ALTITUDE,'geodetic');
Bt=sqrt(BX.^2+BY.^2+BZ.^2);

data_save=[ALTITUDE',Bt];

save('hight_Bvalues.txt','data_save','-ascii');