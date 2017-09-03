clear all;
%反演网格 parameters(均匀网格每个CELL大小相同)
Re=6371.2;    %地球半径或地面台站的位置

Nx=52;
Ny=10;
Nz=40;

%以公里定义;
grid_Latmin=5;
grid_Latmax=109;
d_Lat=(grid_Latmax-grid_Latmin)/Nx;

grid_Lonmin=5;
grid_Lonmax=205;
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


Xgrid_Min = min(Xg);
Xgrid_Max = max(Xg);
Ygrid_Min = min(Yg);
Ygrid_Max = max(Yg);
Zgrid_Min = min(Zg);
Zgrid_Max = max(Zg);

% ray parameters 射线参数
%接收台站 parameters

Re=6371.2;    %地球半径或地面台站的位置

Rec_latmin=0;
Rec_latmax=150;
d_Rec_lat=10;

Rec_lat=Rec_latmin:d_Rec_lat:Rec_latmax;   %沿纬度布了length(rec)个站
Rec_Lon=repmat(70,1,length(Rec_lat));   %经度不变
Rec_Hight=repmat(Re,1,length(Rec_lat));  %地面台站


%卫星 parameters
Rs=6871.2;    %卫星初始高度位置

Sat_latmin=5;
Sat_latmax=79;


Sat_lat=Sat_latmin:1:Sat_latmax;           %方位角信息 75个值走过75个网格 沿着纬度走
Sat_Lon=repmat(70,1,length(Sat_lat));      %卫星经度70公里
Sat_Hight=repmat(Rs,1,length(Sat_lat));    %卫星高度不变


%射线参数 从卫星发射到台站
Nrays = length(Sat_lat)*length(Rec_lat);

% Define 3 dimensional rays
Rays = zeros(Nrays, 2, 3);


x0=reshape(repmat(Sat_lat,length(Rec_lat),1),1,Nrays);
x1=repmat(Rec_lat,1,length(Sat_lat));

y0= reshape(repmat(Sat_Lon,length(Rec_Lon),1),1,Nrays);
y1= repmat(Rec_Lon,1,length(Sat_Lon));

z0= reshape(repmat(Sat_Hight,length(Rec_Hight),1),1,Nrays);
z1= repmat(Rec_Hight,1,length(Sat_Hight));

Rays(:,1,1) = x0;
Rays(:,2,1) = x1;
Rays(:,1,2) = y0;
Rays(:,2,2) = y1;
Rays(:,1,3) = z0;
Rays(:,2,3) = z1;

%% MART coefficient_matrix

[Grid_raysinfo, MART_coefficient_matrix] = MART_coefficient_matrix_3D(Rays,Xg,Yg,Zg);
Grid_intersection_Rays=Grid_raysinfo.Grid_intersection_Rays;
Segment_Lengths_Rays=Grid_raysinfo.Segment_Lengths_Rays;
cells_Rays=Grid_raysinfo.cells_Rays;
MART_coefficient_matrix_3index=Grid_raysinfo.MART_coefficient_matrix_3index;
whos  Grid_intersection_Rays Segment_Lengths_Rays cells_Rays MART_coefficient_matrix_3index MART_coefficient_matrix
save('MART_coefficient_matrix.mat','MART_coefficient_matrix');


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
        axis([Xgrid_Min, Xgrid_Max, Ygrid_Min, Ygrid_Max, Zgrid_Min, Zgrid_Max]);
        grid on
        set(gca, 'XTick', Xg);
        set(gca, 'YTick', Yg);
        set(gca, 'ZTick', Zg);
        hold on
    end
    hold off
end;

%% 用spy画出稀疏矩阵 不同平面的 
%显示系数矩阵
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
%% 显示交点信息
for i = 1:20
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


