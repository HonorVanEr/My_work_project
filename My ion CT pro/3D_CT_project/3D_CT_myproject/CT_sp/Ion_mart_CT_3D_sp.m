%% set up initial parameter
clear all
load('MART_coefficient_matrix_sp.mat');
Nx=52;
Ny=10;
Nz=40;

[Nrays,N_net]=size(MART_coefficient_matrix);
% N_net=2080; %�����ܶ���������
% Nrays=1200; %��������

A=MART_coefficient_matrix*1000;  %���ʵ�λ��
load('Ne2.mat');% ����㾭��ģ��   �߶�xγ��
load('Ne22.mat');%�����Ŷ���ʵ�ʵĵ����ģ��

Ne2=Ne2';
iri_Ne=repmat(Ne2,10,1);
iri_Ne_3d=reshape(iri_Ne,Nx,Ny,Nz);  %����3ά����ʽ
Ne_initial=reshape(iri_Ne_3d,N_net,1);%��ʼֵ ��ʼ�������е����ܶ�����ֲ�
X0=Ne_initial;


Ne22=Ne22';
Ne_true=repmat(Ne22,10,1);
Ne_true_3d=reshape(Ne_true,Nx,Ny,Nz);
Ne_true=reshape(Ne_true,N_net,1);  %ʵ�ʵ����ܶȷֲ�

if 1,
    TEC_R=A*Ne_true;
    Y=TEC_R;
    save('3D_TEC_R_sp.mat','TEC_R');
else
    load('3D_TEC_R_sp.mat');
    Y=TEC_R;
end


if length(Y)==Nrays,
    display('Y��������һ��')
end


%mart(p, A, x0, relax=1., iters=1, tol=0):
%   """Performs multiplicitive algebraic reconstruction (MART) of image`x` given
%     projection data `Y` and a projection matrix `A` which satisfies
%     :math:`\vec{Y} = A\vec{s}`

%     Notes: must assume positivity of image--which makes sense for ionosphere tomography, thus x0 must be > 0
%     x_j^{k+1} = x_j^k * (y_i / (\sum_j a_{ij}x_j^k))^{\gamma\delta_i P_{ij}}



%% begin run
tic

iterations=3;  %��������

% iterations=3;  %��������
%
relax=1;  %�ɳ�����
X= mart(Y,A,X0,relax,iterations);
save('X.mat','X')

toc

%% ��������
Nx=52;
Ny=10;
Nz=40;

Re=6371.2;    %����뾶�����̨վ��λ��

grid_Latmin=14;
grid_Latmax=15.5;
d_Lat=(grid_Latmax-grid_Latmin)/Nx;

grid_Lonmin=70;
grid_Lonmax=80;
d_Lon=(grid_Lonmax-grid_Lonmin)/Ny;

grid_Hightmin=(Re+100);
grid_Hightmax=(Re+500);
d_hight=(grid_Hightmax-grid_Hightmin)/Nz;

%��ʼ������
Lat=grid_Latmin:d_Lat:grid_Latmax;          %����γ�ȷ�Χ -90����γ�� �� 90����γ��
Lon=grid_Lonmin:d_Lon:grid_Lonmax;          %���񾭶ȷ�Χ ����-180 �� ����180
Hight=grid_Hightmin:d_hight:grid_Hightmax;  %����ĸ߶ȷ�Χ ���� 0km ����

Xg=Lat;
Yg=Lon;
Zg=Hight;

%% plot one slice 

% Ne_true_net=reshape(Ne_true,40,52,10);

slice_number=6;
Ne_true_net=Ne_true_3d;% γ�Ⱦ��ȸ߶�

Ne_true_2d=squeeze(Ne_true_net(:,slice_number,:)); %γ�ȸ߶���



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

X_net=reshape(X,Nx,Ny,Nz);   %���ݵõ��ĵ����ܶȷֲ�
X_2d=squeeze(X_net(:,slice_number,:)); %γ�ȸ߶���
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
xslice=Yg(4:6);  % �ؾ���
yslice=[];       %��γ��
zslice=[];       %�ظ߶�

%% PLOT ��ʼ��slice
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


%% PLOT ��ʵ��SLICE 
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


%% PLOT ���ݺ��SLICE 
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




