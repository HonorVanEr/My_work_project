%% set up initial parameter
clear all
load('MART_coefficient_matrix_sp.mat');
Nx=52;
Ny=10;
Nz=40;

[Nrays,N_net]=size(MART_coefficient_matrix);
% N_net=2080; %电子密度网格数量
% Nrays=1200; %射线数量

A=MART_coefficient_matrix*1000;  %国际单位米
load('Ne2.mat');% 电离层经验模型   高度x纬度
load('Ne22.mat');%加上扰动的实际的电离层模型

Ne2=Ne2';
iri_Ne=repmat(Ne2,10,1);
iri_Ne_3d=reshape(iri_Ne,Nx,Ny,Nz);  %化成3维的形式
Ne_initial=reshape(iri_Ne_3d,N_net,1);%初始值 初始估计所有电子密度网格分布
X0=Ne_initial;


Ne22=Ne22';
Ne_true=repmat(Ne22,10,1);
Ne_true_3d=reshape(Ne_true,Nx,Ny,Nz);
Ne_true=reshape(Ne_true,N_net,1);  %实际电子密度分布

if 1,
    TEC_R=A*Ne_true;
    Y=TEC_R;
    save('3D_TEC_R_sp.mat','TEC_R');
else
    load('3D_TEC_R_sp.mat');
    Y=TEC_R;
end


if length(Y)==Nrays,
    display('Y与射线数一样')
end


%mart(p, A, x0, relax=1., iters=1, tol=0):
%   """Performs multiplicitive algebraic reconstruction (MART) of image`x` given
%     projection data `Y` and a projection matrix `A` which satisfies
%     :math:`\vec{Y} = A\vec{s}`

%     Notes: must assume positivity of image--which makes sense for ionosphere tomography, thus x0 must be > 0
%     x_j^{k+1} = x_j^k * (y_i / (\sum_j a_{ij}x_j^k))^{\gamma\delta_i P_{ij}}



%% begin run
tic

iterations=3;  %迭代次数

% iterations=3;  %迭代次数
%
relax=1;  %松弛因子
X= mart(Y,A,X0,relax,iterations);
save('X.mat','X')

toc

%% 定义网格
Nx=52;
Ny=10;
Nz=40;

Re=6371.2;    %地球半径或地面台站的位置

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

Xg=Lat;
Yg=Lon;
Zg=Hight;

%% plot one slice 

% Ne_true_net=reshape(Ne_true,40,52,10);

slice_number=6;
Ne_true_net=Ne_true_3d;% 纬度经度高度

Ne_true_2d=squeeze(Ne_true_net(:,slice_number,:)); %纬度高度面



Ne_true_2d=Ne_true_2d';
v=0:.12e11:1.56e11;
Azi=0:0.0288*82.4891:1.5*82.4891-0.0288*82.4891;
figure
contourf(Azi(9:44)-Azi(9),205:5:400,Ne_true_2d(:,9:44),v);
h=colorbar;
set(get(h,'title'),'string');
axis([0 80 205 400])
xlabel('Azimuth Direcction/km');
ylabel('Altitude/km');
title(['Real Ne' '    slice plane at ' num2str(Lon(slice_number)) '\circ E'],'FontSize',13,'FontWeight','bold');

X_net=reshape(X,Nx,Ny,Nz);   %反演得到的电子密度分布
X_2d=squeeze(X_net(:,slice_number,:)); %纬度高度面
X_2d=X_2d';


v=0:.12e11:1.56e11;
figure
contourf(Azi(9:44)-Azi(9),205:5:400,X_2d(:,9:44),v);
h=colorbar;
set(get(h,'title'),'string');
axis([0 80 205 400])
xlabel('Azimuth Direcction/km');
ylabel('Altitude/km');
title(['CT Inverse Ne' '    slice plane at ' num2str(Lon(slice_number)) '\circ E'],'FontSize',13,'FontWeight','bold');





%% slice choose
xslice=Yg(4:6);  % 沿经度
yslice=[];       %沿纬度
zslice=[];       %沿高度

%% PLOT 初始的slice
figure
[x,y,z] = meshgrid(Yg(1:end-1),Xg(1:end-1),Zg(1:end-1));


h=slice(x,y,z,iri_Ne_3d,xslice,yslice,zslice);

for i=1:length(h)
    h(i).FaceColor = 'interp';
    h(i).EdgeColor = 'none';
    h(i).DiffuseStrength = 0.8;
end
colorbar
caxis([0 1.56e11]);

xlabel('longitudes(\circ)');
ylabel('latitudes(\circ)');
zlabel('Altitude(km)');

title('Initial iri Ne','FontSize',13,'FontWeight','bold');

% h1=gca;
% h1.XLim=minmax(Yg);
% h1.YLim=minmax(Xg);
% h1.ZLim=minmax(Zg);


%% PLOT 真实的SLICE 
figure
[x,y,z] = meshgrid(Yg(1:end-1),Xg(1:end-1),Zg(1:end-1));


h=slice(x,y,z,Ne_true_net,xslice,yslice,zslice);

for i=1:length(h)
    h(i).FaceColor = 'interp';
    h(i).EdgeColor = 'none';
    h(i).DiffuseStrength = 0.8;
end
colorbar
caxis([0 1.56e11]);

xlabel('longitudes(\circ)');
ylabel('latitudes(\circ)');
zlabel('Altitude(km)');

title('Real Ne','FontSize',13,'FontWeight','bold');

% h1=gca;
% h1.XLim=minmax(Yg);
% h1.YLim=minmax(Xg);
% h1.ZLim=minmax(Zg);


%% PLOT 反演后的SLICE 
figure

[x,y,z] = meshgrid(Yg(1:end-1),Xg(1:end-1),Zg(1:end-1));


h=slice(x,y,z,X_net,xslice,yslice,zslice);

for i=1:length(h)
    h(i).FaceColor = 'interp';
    h(i).EdgeColor = 'none';
    h(i).DiffuseStrength = 0.8;
end
colorbar
caxis([0 1.56e11]);

xlabel('longitudes(\circ)');
ylabel('latitudes(\circ)');
zlabel('Altitude(km)');

title('3D CT Inverse Ne','FontSize',13,'FontWeight','bold');

% h1=gca;
% h1.XLim=minmax(Yg);
% h1.YLim=minmax(Xg);
% h1.ZLim=minmax(Zg);




