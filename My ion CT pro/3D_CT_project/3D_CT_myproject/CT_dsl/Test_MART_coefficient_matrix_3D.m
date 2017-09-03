%% 测试参数  测试求直角坐标系下的系数矩阵
clear all
% define grid
Xg = 10*[-2 1 3 5 6 8 10 13];
Yg = 10*[-1 1 3 5 6 8 10];
Zg = 10*[ 0 1 3 5 6 8];
Xgrid_Min = min(Xg);
Xgrid_Max = max(Xg);
Ygrid_Min = min(Yg);
Ygrid_Max = max(Yg);
Zgrid_Min = min(Zg);
Zgrid_Max = max(Zg);

% define rays
Nrays = 4;
Ray_Ind = 1:Nrays;

% Define 3 dimensional rays
x0 = Xgrid_Min + 15*Ray_Ind;
x1 = Xgrid_Max - 8*Ray_Ind;
y0 = Ygrid_Min + 10 + Ray_Ind;
y1 = Ygrid_Max - 20 - Ray_Ind;
z0 = Zgrid_Min + 20 + Ray_Ind;
z1 = Zgrid_Max - 10 - Ray_Ind;
%x1 = x0;
%y1 = y0;
%z1 = z0;

% 3 dimensional example
Rays = zeros(Nrays, 2, 3);
Rays(:,1,1) = x0;
Rays(:,2,1) = x1;
Rays(:,1,2) = y0;
Rays(:,2,2) = y1;
Rays(:,1,3) = z0;
Rays(:,2,3) = z1;

% [Segment_Lengths, Cells, Valid_Intersections, Intersections] = find_ray_grid_intersections(Rays, Xg, Yg, Zg);

%% MART coefficient_matrix
[Grid_raysinfo, MART_coefficient_matrix] = MART_coefficient_matrix_3D(Rays,Xg,Yg,Zg);
Grid_intersection_Rays=Grid_raysinfo.Grid_intersection_Rays;
Segment_Lengths_Rays=Grid_raysinfo.Segment_Lengths_Rays;
cells_Rays=Grid_raysinfo.cells_Rays;
MART_coefficient_matrix_3index=Grid_raysinfo.MART_coefficient_matrix_3index;

%% plot Rays Intersections
plot_display=1;

if plot_display,
    figure
    for i = 1:Nrays
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


