time = datenum([2007 7 17 6 30 0]);
lat_start = 42.17; % Geodetic latitudes in degrees.
lon_start = 128.0; % Geodetic longitudes in degrees.
Re=6371.2; %km
ALTITUDE =205:5:400; % Altitude in km.


cd /Users/yangjian/Documents/MYgithub/Matlab_Geospace_Analysis_Package/gap/gap_library/other_toolboxes/igrf_2012
load igrfcoefs.mat
[BX, BY, BZ]=igrf(time, lat_start, lon_start, ALTITUDE,'geodetic');
Bt=sqrt(BX.^2+BY.^2+BZ.^2);

data_save=[ALTITUDE',Bt];


cd /Users/yangjian/Desktop
save('hight_Bvalues_igrf2012.txt','data_save','-ascii');

if 1,
    
    % use geocentric
    R=Re+ALTITUDE;
    [BX, BY, BZ]=igrf(time, lat_start, lon_start, R,'geocentric');
    Bt=sqrt(BX.^2+BY.^2+BZ.^2);
    
    data_save=[ALTITUDE',Bt];
    
    save('hight_Bvalues_igrf2012_geocentric.txt','data_save','-ascii');
    
end


cd /Users/yangjian/Documents/MYgithub/Matlab_Geospace_Analysis_Package/gap/gap_library/other_toolboxes/igrf_2012
