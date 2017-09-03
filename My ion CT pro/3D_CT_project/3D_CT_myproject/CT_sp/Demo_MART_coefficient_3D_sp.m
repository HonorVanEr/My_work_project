%% ���巴���������� define grid
clear all;

%geodetic
%�������� parameters(��������ÿ��CELL��С��ͬ)
Re=6371.2;    %����뾶�����̨վ��λ��

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

%��ʼ������
Lat=grid_Latmin:d_Lat:grid_Latmax;          %����γ�ȷ�Χ -90����γ�� �� 90����γ��
Lon=grid_Lonmin:d_Lon:grid_Lonmax;          %���񾭶ȷ�Χ ����-180 �� ����180
Hight=grid_Hightmin:d_hight:grid_Hightmax;  %����ĸ߶ȷ�Χ ���� 0km ����

Xg=Lat;
Yg=Lon;
Zg=Hight;

%% ray parameters ���߲��� ��γ������ϵ
%����̨վ parameters

Re=6371.2;    %����뾶�����̨վ��λ�� km

%���沿�����վ ��γ��
Rec_latmin=14.0;
Rec_latmax=15.9;
d_Rec_lat=0.2;

Rec_lat=Rec_latmin:d_Rec_lat:Rec_latmax;   %��γ�Ȳ���length(Rec_lat)��վ
% Rec_Lon=repmat(70,1,length(Rec_lat));   %���Ȳ���
lon_set=75;
Rec_pos_lat=[Rec_lat;repmat(lon_set,1,length(Rec_lat));repmat(Re,1,length(Rec_lat))];

%���沿�����վ �ؾ���
Rec_lonmin=70;
Rec_lonmax=80;
d_Rec_lon=2;

Rec_Lon=Rec_lonmin:d_Rec_lon:Rec_lonmax;   %���Ȳ�length(Rec_Lon)��̨վ
lat_set=14.5;

Rec_pos_lon=[repmat(lat_set,1,length(Rec_Lon));Rec_Lon;repmat(Re,1,length(Rec_Lon))];

%���沿�����վ �ظ߶�
Rec_Hight=repmat(Re,1,length(Rec_lat));  %���Ÿ߶Ȳ�̨վ ����γ�ȶ�Ϊ��ֵ
Rec_pos_hight=[repmat(0,1,length(Rec_Hight));repmat(0,1,length(Rec_Hight));Rec_Hight];
Rec_pos_hight=[];

%���е������̨վλ�� ѡ���һ�޶���̨վ
Rec_pos=[Rec_pos_lat Rec_pos_lon Rec_pos_hight];
[Rec_pos,~,~] = unique(Rec_pos','rows','stable');



%���� parameters ����һ�������еĸ߶�Ϊ��ֵ�����������ؾ�����
Rs=6871.2;    %���ǳ�ʼ�߶�λ��

G=6.67*10^(-12); %������������
Mz=5.976*10^24;
T=2*pi*(6871200^1.5)/(G*Mz)^0.5;%�����˶����ڿ����ն���
f=1/T*360;               %γ����Ϣ

Sat_latmin=14.0;
Sat_latmax=15.5;

Sat_lat=Sat_latmin:f:Sat_latmax;           %��λ����Ϣ 75��ֵ�߹�75������ ����γ����
Sat_Lon=repmat(75,1,length(Sat_lat));      %���Ǿ��Ȳ��� E 75
Sat_Hight=repmat(Rs,1,length(Sat_lat));    %���Ǹ߶Ȳ���

Sat_pos=[Sat_lat;Sat_Lon;Sat_Hight];

[Sat_pos,~,~] = unique(Sat_pos','rows','stable');



%���߲��� �����Ƿ��䵽̨վ
Nrays = length(Sat_pos(:,1))*length(Rec_pos(:,1));

% Define 3 dimensional rays
Rays = zeros(Nrays, 2, 3);



temp1=reshape(repmat(Rec_pos(:,1)',1,length(Sat_pos(:,1))),1,Nrays);
temp2=reshape(repmat(Rec_pos(:,2)',1,length(Sat_pos(:,1))),1,Nrays);
temp3=reshape(repmat(Rec_pos(:,3)',1,length(Sat_pos(:,1))),1,Nrays);
rays_start=[temp1' temp2' temp3'];

temp1=reshape(repmat(Sat_pos(:,1)',length(Rec_pos(:,1)),1),1,Nrays);
temp2=reshape(repmat(Sat_pos(:,2)',length(Rec_pos(:,1)),1),1,Nrays);
temp3=reshape(repmat(Sat_pos(:,3)',length(Rec_pos(:,1)),1),1,Nrays);
rays_end=[temp1' temp2' temp3'];

ecef=0;
if ecef,
    %Convert geodetic coordinates to ECEF coordinates
    [x0,y0,z0]=geod2ecef(rays_start(:,1), rays_start(:,2), (rays_start(:,3)-Re)*1000);
    
    [x1,y1,z1]=geod2ecef(rays_end(:,1), rays_end(:,2), (rays_end(:,3)-Re)*1000);
    
    %Convert geodetic spherical coordinates to  cart coordinates
else
    [x0,y0,z0]=sph2cart(degtorad(rays_start(:,2)), degtorad(rays_start(:,1)), rays_start(:,3)*1000);
    
    [x1,y1,z1]=sph2cart(degtorad(rays_end(:,2)), degtorad(rays_end(:,1)), rays_end(:,3)*1000);
    
end
%����ת��
% γ�ȷ�Χ -90����γ�� �� 90����γ��
% ���ȷ�Χ ����-180 �� ����180
% �߶ȷ�Χ ���� 0km ����



Rays(:,1,1) = x0;
Rays(:,2,1) = x1;
Rays(:,1,2) = y0;
Rays(:,2,2) = y1;
Rays(:,1,3) = z0;
Rays(:,2,3) = z1;

% save('Rays.mat','Rays');

%% MART coefficient_matrix

[Grid_raysinfo, MART_coefficient_matrix] = MART_coefficient_matrix_3D_sp(Rays,Xg,Yg,Zg);
Grid_intersection_Rays=Grid_raysinfo.Grid_intersection_Rays;
Grid_intersection_sp=Grid_raysinfo.Grid_intersection_sp;

Segment_Lengths_Rays=Grid_raysinfo.Segment_Lengths_Rays;
cells_Rays=Grid_raysinfo.cells_Rays;
MART_coefficient_matrix_3index=Grid_raysinfo.MART_coefficient_matrix_3index;
whos  Grid_intersection_Rays Segment_Lengths_Rays cells_Rays MART_coefficient_matrix_3index MART_coefficient_matrix
save('MART_coefficient_matrix_sp.mat','MART_coefficient_matrix');

%% plot Rays Intersections
plot_display=1;
dNrays=100;
if plot_display,
    figure
    for i = 1:dNrays:Nrays
        x0=Rays(i,1,1);
        x1=Rays(i,2,1);
        y0=Rays(i,1,2);
        y1=Rays(i,2,2);
        z0=Rays(i,1,3);
        z1=Rays(i,2,3);
        Rays_x = [x0;x1];
        Rays_y = [y0;y1];
        Rays_z = [z0;z1];
        Xi = Grid_intersection_Rays{i,:}(:,1);
        Yi = Grid_intersection_Rays{i,:}(:,2);
        Zi = Grid_intersection_Rays{i,:}(:,3);
        
        plot3(Rays_x, Rays_y, Rays_z,Xi,Yi,Zi,'*');
        title('3D Ray Intersections with Irregular 3-D Rectangular Grid');
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        %axis([Xgrid_Min, Xgrid_Max, Ygrid_Min, Ygrid_Max, Zgrid_Min, Zgrid_Max]);
        grid on
        %set(gca, 'XTick', Xg);
        %set(gca, 'YTick', Yg);
        %set(gca, 'ZTick', Zg);
        hold on
    end
    hold off
end;


%��ʾϵ������
plot_spy=1;
trace_number=30;
if plot_spy,
    mart_temp=reshape(MART_coefficient_matrix_3index(trace_number,:,:,:),(length(Xg)-1),(length(Yg)-1),(length(Zg)-1));
    
    figure(2);
    title('YZ plane Intercepts');
    hold on
    for k=1:length(Xg)-1
        spy(reshape(mart_temp(k,:,:),(length(Yg)-1),(length(Zg)-1)));
    end
    hold off
    
    figure(3);
    title('XZ plane Intercepts');
    hold on
    for k=1:length(Yg)-1
        spy(reshape(mart_temp(:,k,:),(length(Xg)-1),(length(Zg)-1)));
    end
    hold off
    
    figure(4);
    title('XY plane Intercepts');
    hold on
    for k=1:length(Zg)-1
        spy(reshape(mart_temp(:,:,k),(length(Xg)-1),(length(Yg)-1)));
    end
    hold off
end
%% ��ʾ������Ϣ
for i = 1:10
    if size(Segment_Lengths_Rays{i,:},1)>0,
        fprintf('Ray %d Segment Lengths\n', i);
        s = find(Segment_Lengths_Rays{i,:} > 0);
        
        for j = 1:length(s)
            fprintf('X cell = %2d, Y cell = %2d, Z cell = %2d, Length = %15.8f\n',...
                cells_Rays{i,:}(j,1),cells_Rays{i,:}(j,2),cells_Rays{i,:}(j,3), Segment_Lengths_Rays{i,:}(j));
        end;
        
        fprintf('Ray %d Grid Intercepts\n', i);
        for j = 1:length(Grid_intersection_Rays{i,:}(:,1))
            fprintf('Xi = %5.1f, Yi = %5.1f, Zi = %5.1f\n',...
                Grid_intersection_Rays{i,:}(j,1), Grid_intersection_Rays{i,:}(j,2),Grid_intersection_Rays{i,:}(j,3));
        end;
    else
        fprintf('Ray %d NO Segment Lengths\n', i);
        
    end
    
end


