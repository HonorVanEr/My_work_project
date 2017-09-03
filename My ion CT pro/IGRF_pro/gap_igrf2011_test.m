%% Test gap_igrf2011
cd /Users/yangjian/Documents/MYgithub/Matlab_Geospace_Analysis_Package/gap/gap_library/gap_igrf2011/
gh_coeffs = gap_prepare_igrf();


fyears=2007.5425;
itype  = 2 ;  %geocentric (sphere)

alt=6371.2+205; %km
nlat=42.17;%deg
elong=128.0;%deg

%%  IGRF magnetic  geocentric
[B]=crino_igrf11syn(gh_coeffs,fyears,itype, alt,nlat,elong)
B_t=sqrt(B(:,1)^2+B(:,2)^2+B(:,3)^2)


%% Gives the IGRF magnetic field model in specified coordinates
gh_coeffs = gap_prepare_igrf();
times=datenum(2007,07,17);
r=6371.2+205; %km
lat=42.17;%deg
long=128.0;%deg


xyz=[long, lat, r];
xyz_frame='NEC';
end_frame='GEO';
itype=2;  %geocentric (sphere)

[B_end_frame] = gap_igrf(gh_coeffs, times, xyz, xyz_frame, end_frame, itype)

B_t2=sqrt(B_end_frame(:,1)^2+B_end_frame(:,2)^2+B_end_frame(:,3)^2)










