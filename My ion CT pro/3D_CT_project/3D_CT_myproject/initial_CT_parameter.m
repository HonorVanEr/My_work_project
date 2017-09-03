%% 定义反演网格区域 define grid
clear all

%反演网格 parameters(均匀网格每个CELL大小相同)
Re=6371.2;    %地球半径或地面台站的位置

Nx=52;
Ny=10;
Nz=40;

grid_Latmin=14;
grid_Latmax=15.5;
d_Lat=(grid_Latmax-grid_Latmin)/Nx;

grid_Lonmin=70;
grid_Lonmax=80;
d_Lon=(grid_Lonmax-grid_Lonmin)/Ny;

grid_Hightmin=(Re+100);
grid_Hightmax=(Re+500);
d_hight=(grid_Hightmax-grid_Hightmin)/Nz;

%初始化网格
Lat=grid_Latmin:d_Lat:grid_Latmax;          %网格纬度范围 -90（南纬） 到 90（北纬）
Lon=grid_Lonmin:d_Lon:grid_Lonmax;          %网格经度范围 西经-180 到 东经180
Hight=grid_Hightmin:d_hight:grid_Hightmax;  %网格的高度范围 地面 0km 以上

% define grid
Lon_temp=interp1(0:Ny,Lon,0:Ny/Nx:Ny);
Hight_temp=interp1(0:Nz,Hight,0:Nz/Nx:Nz);



[Xg,Yg,Zg]=sph2cart(degtorad(Lon_temp),degtorad(Lat),Hight_temp);
Xgrid_Min = min(Xg);
Xgrid_Max = max(Xg);
Ygrid_Min = min(Yg);
Ygrid_Max = max(Yg);
Zgrid_Min = min(Zg);
Zgrid_Max = max(Zg);

Xg_near=Xgrid_Min:(Xgrid_Max-Xgrid_Min)/Nx:Xgrid_Max;
Yg_near=Ygrid_Min:(Ygrid_Max-Ygrid_Min)/Ny:Ygrid_Max;
Zg_near=Zgrid_Min:(Zgrid_Max-Zgrid_Min)/Nz:Zgrid_Max;




%% ray parameters 射线参数 经纬度坐标系
%接收台站 parameters

Re=6371.2;    %地球半径或地面台站的位置

Rec_latmin=14.6;
Rec_latmax=14.9;
d_Rec_lat=0.02;

Rec_lat=Rec_latmin:d_Rec_lat:Rec_latmax;   %沿纬度布了length(rec)个站
Rec_Lon=repmat(70,1,length(Rec_lat));   %经度不变
Rec_Hight=repmat(Re,1,length(Rec_lat));  %地面台站


%卫星 parameters
Rs=6871.2;    %卫星初始高度位置

G=6.67*10^(-12); %万有引力常量
Mz=5.976*10^24;
T=2*pi*(6871200^1.5)/(G*Mz)^0.5;%卫星运动周期开普勒定律
f=1/T*360;               %纬度信息

Sat_latmin=14;
Sat_latmax=15.5;

Sat_lat=Sat_latmin:f:Sat_latmax;           %方位角信息 75个值走过75个网格 沿着纬度走
Sat_Lon=repmat(70,1,length(Sat_lat));      %卫星经度不变
Sat_Hight=repmat(Rs,1,length(Sat_lat));    %卫星高度不变


%射线参数 从卫星发射到台站
Nrays = length(Sat_lat)*length(Rec_lat);

% Define 3 dimensional rays
Rays = zeros(Nrays, 2, 3);

x0_sph=reshape(repmat(Sat_lat,length(Rec_lat),1),1,Nrays);
x1_sph=repmat(Rec_lat,1,length(Sat_lat));

y0_sph = reshape(repmat(Sat_Lon,length(Rec_Lon),1),1,Nrays);
y1_sph = repmat(Rec_Lon,1,length(Sat_Lon));

z0_sph = reshape(repmat(Sat_Hight,length(Rec_Hight),1),1,Nrays);
z1_sph = repmat(Rec_Hight,1,length(Sat_Hight));


%坐标转化
% 纬度范围 -90（南纬） 到 90（北纬）
% 经度范围 西经-180 到 东经180
% 高度范围 地面 0km 以上
[Rec_x,Rec_y,Rec_z]=sph2cart(degtorad(Rec_Lon),degtorad(Rec_lat),Rec_Hight);
[Sat_x,Sat_y,Sat_z]=sph2cart(degtorad(Sat_Lon),degtorad(Sat_lat),Sat_Hight);


[x0,y0,z0]=sph2cart(degtorad(y0_sph),degtorad(x0_sph),z0_sph);
[x1,y1,z1]=sph2cart(degtorad(y1_sph),degtorad(x1_sph),z1_sph);

Rays(:,1,1) = x0;
Rays(:,2,1) = x1;
Rays(:,1,2) = y0;
Rays(:,2,2) = y1;
Rays(:,1,3) = z0;
Rays(:,2,3) = z1;

% save('Rays.mat','Rays');

% [Segment_Lengths, Cells, Valid_Intersections, Intersections] = find_ray_grid_intersections(Rays, Xg, Yg, Zg);

%% MART coefficient_matrix


MART_coefficient_matrix=zeros(length(Rays),length(Lat),length(Lon),length(Hight));






