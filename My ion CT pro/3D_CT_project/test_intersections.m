% Find locations where a multiple rays intersect irregular N-dimensional grid lines
% Base this on a parametric form of the line equation
clear
close all
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
Nrays = 6;
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

[Segment_Lengths, Cells, Valid_Intersections, Intersections] = find_ray_grid_intersections(Rays, Xg, Yg, Zg);

Rays_x = [x0;x1];
Rays_y = [y0;y1];
Rays_z = [z0;z1];

Xi = Intersections(:,:, 1);
Yi = Intersections(:,:, 2);
Zi = Intersections(:,:, 3);

figure
plot3(Rays_x, Rays_y, Rays_z, Xi(Valid_Intersections), Yi(Valid_Intersections), Zi(Valid_Intersections), '*');
title('3D Ray Intersections with Irregular 3-D Rectangular Grid');
xlabel('X');
ylabel('Y');
zlabel('Z');
axis([Xgrid_Min, Xgrid_Max, Ygrid_Min, Ygrid_Max, Zgrid_Min, Zgrid_Max]);
grid on
set(gca, 'XTick', Xg);
set(gca, 'YTick', Yg);
set(gca, 'ZTick', Zg);

for i = 1:Nrays
    fprintf('Ray %d Segment Lengths\n', i);
    s = find(Segment_Lengths(i,:) > 0);
    for j = 1:length(s)
        k = s(j);
        fprintf('X cell = %2d, Y cell = %2d, Z cell = %2d, Length = %5.1f\n', Cells(i, k, 1), Cells(i, k, 2), Cells(i, k, 3), Segment_Lengths(i,k));
    end;
    fprintf('Ray %d Grid Intercepts\n', i);
        s = find(Valid_Intersections(i,:) > 0);
    for j = 1:length(s)
        k = s(j);
        fprintf('Xi = %5.1f, Yi = %5.1f, Zi = %5.1f\n', Xi(i,k), Yi(i,k), Zi(i,k));
    end;

end

% 2-dimensional example
Rays = zeros(Nrays, 2, 2);
Rays(:,1,1) = x0;
Rays(:,2,1) = x1;
Rays(:,1,2) = y0;
Rays(:,2,2) = y1;

[Segment_Lengths, Cells, Valid_Intersections, Intersections] = find_ray_grid_intersections(Rays, Xg, Yg);

Xi = Intersections(:,:, 1);
Yi = Intersections(:,:, 2);

figure
plot(Rays_x, Rays_y, Xi(Valid_Intersections), Yi(Valid_Intersections), '*');
title('2D Ray Intersections with Irregular 2-D Rectangular Grid');
xlabel('X');
ylabel('Y');
axis([Xgrid_Min, Xgrid_Max, Ygrid_Min, Ygrid_Max]);
grid on
set(gca, 'XTick', Xg);
set(gca, 'YTick', Yg);

for i = 1:Nrays
    fprintf('Ray %d Segment Lengths\n', i);
    s = find(Segment_Lengths(i,:) > 0);
    for j = 1:length(s)
        k = s(j);
        fprintf('X cell = %2d, Y cell = %2d, Length = %5.1f\n', Cells(i, k, 1), Cells(i, k, 2), Segment_Lengths(i,k));
    end;
    fprintf('Ray %d Grid Intercepts\n', i);
        s = find(Valid_Intersections(i,:) > 0);
    for j = 1:length(s)
        k = s(j);
        fprintf('Xi = %5.1f, Yi = %5.1f\n', Xi(i,k), Yi(i,k));
    end;

end