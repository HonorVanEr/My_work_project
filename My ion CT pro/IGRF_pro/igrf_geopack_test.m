%% 用的都是地心坐标系 不太精确


addpath('/Users/yangjian/Documents/MYgithub/My_matlab_pro/jian_matlab_lib/GEOPACK');
clc;
clear all;

GEOPACK_RECALC(2007,198,6,30,0);

Re=6371.2; %km
hight=205:5:400;
R=(Re+hight)/Re;
THETA=zeros(1,length(R));
THETA(:)=42.17; % Geodetic latitudes in degrees.

PHI=zeros(1,length(R));
PHI(:)= 128; % Geodetic longitudes in degrees.


%角度转弧度制
THETA=THETA*pi/180.0;
PHI=PHI*pi/180.0;

%%

for i=1:length(R)
[BR(i),BTHETA(i),BPHI(i)] = GEOPACK_IGRF_GEO(R(i),THETA(i),PHI(i));
end

Bt=sqrt(BR.^2+BTHETA.^2+BPHI.^2);

data_save=[hight',Bt'];

save('Geopack_hight_Bvalues.txt','data_save','-ascii');