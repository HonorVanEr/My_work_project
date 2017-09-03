%% 测试网格射线参数
clear all;
Xg=[14.5:0.1:15.5];
Yg=[70:1:80];
Hight=[100:20:500];%距离地面的高度 km
Re=6371.2;
Zg=Re+Hight;
rays_start=[15.0,72,6371.2;15.1,75,6371.2];
rays_end=[15.1,73,6371.2+550;15.2,78,6371.2+550];

[x0,y0,z0]=sph2cart(degtorad(rays_start(:,2)), degtorad(rays_start(:,1)),rays_start(:,3)*1000);
[x1,y1,z1]=sph2cart(degtorad(rays_end(:,2)), degtorad(rays_end(:,1)), rays_end(:,3)*1000);
Nrays=length(rays_start(:,3));
Rays = zeros(Nrays, 2, 3);
Rays(:,1,1) = x0;
Rays(:,2,1) = x1;
Rays(:,1,2) = y0;
Rays(:,2,2) = y1;
Rays(:,1,3) = z0;
Rays(:,2,3) = z1;


%% 计算系数矩阵
[Grid_raysinfo, MART_coefficient_matrix]=MART_coefficient_matrix_3D_sp(Rays,Xg,Yg,Zg);
Grid_intersection_Rays=Grid_raysinfo.Grid_intersection_Rays;
Grid_intersection_sp=Grid_raysinfo.Grid_intersection_sp;

Segment_Lengths_Rays=Grid_raysinfo.Segment_Lengths_Rays;
cells_Rays=Grid_raysinfo.cells_Rays;
MART_coefficient_matrix_3index=Grid_raysinfo.MART_coefficient_matrix_3index;
whos  Grid_intersection_Rays Segment_Lengths_Rays cells_Rays MART_coefficient_matrix_3index MART_coefficient_matrix

%% plot Rays Intersections
plot_display=1;
dNrays=1;
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


%% 显示系数矩阵
plot_spy=1;
trace_number=1;
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
for i = 1:Nrays
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














