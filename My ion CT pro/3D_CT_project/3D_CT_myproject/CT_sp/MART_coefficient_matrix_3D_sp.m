function [Grid_raysinfo, MART_coefficient_matrix] = MART_coefficient_matrix_3D_sp(Rays,Xg,Yg,Zg)

%计算直角坐标系下多条射线与空间三维网格的交点以及每个网格的截距
%Calculate the length of ray segments within an n-dimensional irregular grid
%   and also intersection locations for multiple rays
%
%
% INPUTS:
%   Rays contains all rays that may intersect the grid
%   Rays is an No_Rays by 2 by No_Dims 3-dimensional array
%   No_Rays is the number of rays
%   No_Dims is the number of dimensions
%   Rays(i, 1, :) contains the starting coordinates for ray i
%   Rays(i, 2, :) contains the ending   coordinates for ray i
%   varargin contains No_Dims vectors. Each vector contains the grid
%   locations for that coordinate.
%   Rays Nx2X3     N 射线的条数
%   Rays(i, 1, :)  第i条射线的起点坐标 台站位置
%   Rays(i, 2, :)  第i条射线的终点坐标 卫星位置
% 网格参数
%   Xg    lat 纬度格点  纬度范围 -90（南纬） 到 90（北纬）deg
%   Yg    lon 经度格点  经度范围 西经-180 到 东经180 deg
%   Zg    高度 高度格点  网格的高度范围 离地心的距离 km
%
% OUTPUTS:
%      Grid_raysinfo 结构体
%        Grid_i ntersection_Rays 交点坐标直角坐标系
%        Grid_intersection_sp   交点坐标球坐标系
%        Segment_Lengths_Rays 射线截取的线段截距
%        cells_Rays    截距的三维下标
%        MART_coefficient_matrix_3index 系数矩阵 N*Nx*Ny*Nz



%      MART_coefficient_matrix 系数矩阵   Nx(Nx*Ny*Nz) 单位km
%
% Examples:
% Xg=[14.5:0.1:15.5];
% Yg=[70:1:80];
% Hight=[100:20:500];%距离地面的高度 km
% Re=6371.2;
% Zg=Re+Hight;
% rays_start=[15.0,72,6371.2;15.1,75,6371.2];
% rays_end=[15.1,73,6371.2+550;15.2,78,6371.2+550];
%
% [x0,y0,z0]=sph2cart(degtorad(rays_start(:,2)), degtorad(rays_start(:,1)),rays_start(:,3)*1000);
% [x1,y1,z1]=sph2cart(degtorad(rays_end(:,2)), degtorad(rays_end(:,1)), rays_end(:,3)*1000);
% Nrays=length(rays_start(:,3));
% Rays = zeros(Nrays, 2, 3);
% Rays(:,1,1) = x0;
% Rays(:,2,1) = x1;
% Rays(:,1,2) = y0;
% Rays(:,2,2) = y1;
% Rays(:,1,3) = z0;
% Rays(:,2,3) = z1;



% [Grid_raysinfo, MART_coefficient_matrix]=MART_coefficient_matrix_3D_sp(Rays,Xg,Yg,Zg);
%

% Ion_MART coefficient_matrix
%
Re=6371.2;    %地球半径或地面台站的位置


Nrays=size(Rays,1);
N_net=(length(Xg)-1)*(length(Yg)-1)*(length(Zg)-1);
N_gridx=length(Xg);
N_gridy=length(Yg);
N_gridz=length(Zg);
Xgrid_Min = min(Xg);
Xgrid_Max = max(Xg);
Ygrid_Min = min(Yg);
Ygrid_Max = max(Yg);
Zgrid_Min = min(Zg);
Zgrid_Max = max(Zg);

d_Lat=Xg(2)-Xg(1);
d_Lon=Yg(2)-Yg(1);
d_hight=Zg(2)-Zg(1);
Grid_intersection_Rays=cell(Nrays,1);

MART_coefficient_matrix_3index=zeros(Nrays,(length(Xg)-1),(length(Yg)-1),(length(Zg)-1));
MART_coefficient_matrix=zeros(Nrays,N_net);

for i=1:Nrays;
    x0=Rays(i,1,1);
    x1=Rays(i,2,1);
    y0=Rays(i,1,2);
    y1=Rays(i,2,2);
    z0=Rays(i,1,3);
    z1=Rays(i,2,3);
    
    
    line_vet=[x1-x0 y1-y0 z1-z0];
    line_point=[x0 y0 z0];
    
    SC_point=[x1 y1 z1];
    GS_point=[x0 y0 z0];
    
    height=(Zg-Re)*1000; %m
    %高度平面交点
    [x_height,y_height,z_height]=Ion_line_plane_intersect(SC_point,GS_point,'height',height);
    xy_intersection=[x_height;y_height;z_height]';
    if size(xy_intersection,2) ~= 3,
        xy_intersection=reshape(xy_intersection,length(xy_intersection)/3,3);
    end
    
    %检验
    %     [az,el,r] = cart2sph(xy_intersection(:,1),xy_intersection(:,2),xy_intersection(:,3));
    %     lat_test=rad2deg(el); lon_test=rad2deg(az);
    
    if isempty(xy_intersection),xy_intersection=[];end
    
    %经度平面交点
    [x_lon,y_lon,z_lon]=Ion_line_plane_intersect(SC_point,GS_point,'lon',Yg); %#ok<*ASGLU>
    xz_intersection=[x_lon;y_lon;z_lon]';
    if size(xz_intersection,2) ~= 3,
        xz_intersection=reshape(xz_intersection,length(xz_intersection)/3,3);
    end
    
    if isempty(xz_intersection),xz_intersection=[];end
    
    %纬度平面交点
    [x_lat,y_lat,z_lat]=Ion_line_plane_intersect(SC_point,GS_point,'lat',Xg);
    yz_intersection=[x_lat;y_lat;z_lat]';
    if size(yz_intersection,2) ~= 3,
        yz_intersection=reshape(yz_intersection,length(yz_intersection)/3,3);
    end
    
    
    if isempty(yz_intersection),yz_intersection=[];end
    
    
    %所有交点坐标排序理论上按一个维度排列即可 排除相同交点
    Grid_intersection=[xy_intersection;xz_intersection;yz_intersection];
    [Grid_intersection, ~, ~] = unique(Grid_intersection,'rows','stable');    %自动升序排序了？
    
    %剔除nan值
    Grid_intersection(any(isnan(Grid_intersection),2),:)=[];
    
    if isempty(Grid_intersection),
        continue
    end
    
    
    %沿X方向排序
    if line_vet(1)>0,
        temp=sortrows(Grid_intersection,1);
    else
        temp=sortrows(Grid_intersection,-1);
    end
    %沿Y方向排序
    if line_vet(2)>0,
        temp2=sortrows(temp,2);
    else
        temp2=sortrows(temp,-2);
    end
    %沿Z方向排序
    if line_vet(3)>0,
        temp3=sortrows(temp2,3);
    else
        temp3=sortrows(temp2,-3);
    end
    
    %break
    Grid_intersection=temp3;
    clear temp temp2 temp3
    
    
    %找到位于射线上的点
    x_ind=find(Grid_intersection(:,1)>=min([x0 x1]) & Grid_intersection(:,1)<=max([x0 x1]));
    y_ind=find(Grid_intersection(:,2)>=min([y0 y1]) & Grid_intersection(:,2)<=max([y0 y1]));
    z_ind=find(Grid_intersection(:,3)>=min([z0 z1]) & Grid_intersection(:,3)<=max([z0 z1]));
    Grid_intersection=Grid_intersection(intersect(intersect(x_ind,y_ind),z_ind),:);
    
    %找到位于网格内的点 两种直角坐标系网格和 经纬高坐标系网格
    [az,el,r] = cart2sph(Grid_intersection(:,1),Grid_intersection(:,2),Grid_intersection(:,3));
    
    if verLessThan('matlab','9.1.0'),
        latitude=radtodeg(el); longitude=radtodeg(az); altitude=r;
    else
        latitude=rad2deg(el); longitude=rad2deg(az); altitude=r;
    end
    
    x_ind=find(latitude>=Xgrid_Min & latitude<=Xgrid_Max);
    y_ind=find(longitude>=Ygrid_Min & longitude<=Ygrid_Max);
    z_ind=find(altitude/1000.0>=Zgrid_Min & altitude/1000.0<=Zgrid_Max);
    Grid_intersection=Grid_intersection(intersect(intersect(x_ind,y_ind),z_ind),:);
    
    
    %距离非常近的认为是一个点
    tol=1.0e-2;
    
    Grid_intersection_Rays{i,1}=Grid_intersection;
    if (size(Grid_intersection,1)>1),
        dRi = diff(Grid_intersection, 1);
        Segment_Lengths = sqrt(sum(dRi.^2, 2));
        
        if find(Segment_Lengths<=tol)
            Grid_intersection(find(Segment_Lengths<=tol)+1,:)=[];
            Segment_Lengths(find(Segment_Lengths<=tol))=[];
        end
        Segment_Lengths_Rays{i,1}=Segment_Lengths; %#ok<SAGROW>
    else
        Segment_Lengths_Rays{i,1}={};
        continue
    end
    
    
    %[latitude, longitude, altitude] = ecef2geod(Grid_intersection(:,1),Grid_intersection(:,2),Grid_intersection(:,3));
    
    
    [az,el,r] = cart2sph(Grid_intersection(:,1),Grid_intersection(:,2),Grid_intersection(:,3));
    if verLessThan('matlab','9.1.0'),
        latitude=radtodeg(el); longitude=radtodeg(az); altitude=r;
    else
        latitude=rad2deg(el); longitude=rad2deg(az); altitude=r;
    end
    clear az el r;
    Grid_intersection_sp=[latitude longitude altitude/1.0e3];
    
    
    %交点纬度经度不是严格的升序或降序
    for j=1:length(Segment_Lengths);
        %判断两交点中点位置所在cell 两种方法
        
        if 1,
            x_cell_index=floor(((latitude(j)+latitude(j+1))/2-Xg(1))/d_Lat)+1;
            y_cell_index=floor(((longitude(j)+longitude(j+1))/2-Yg(1))/d_Lon)+1;
            z_cell_index=floor(((altitude(j)+altitude(j+1))/1000.0/2-Zg(1))/d_hight)+1;
            
            
        else
            x_cell_index=find(Xg<=round((latitude(j)+latitude(j+1))/2,5));
            if isempty(x_cell_index),
                x_cell_index=1;
            end
            x_cell_index=x_cell_index(end);
            if x_cell_index==N_gridx,
                x_cell_index=x_cell_index-1;
            end
            
            y_cell_index=find(Yg<=round((longitude(j)+longitude(j+1))/2,5));
            y_cell_index=y_cell_index(end);
            if y_cell_index==N_gridy,
                y_cell_index=y_cell_index-1;
            end
            
            %经度不够用四舍五入
            z_cell_index=find(Zg<=round((altitude(j)+altitude(j+1))/2000.0,5));
            if isempty(z_cell_index),
                z_cell_index=1;
            end
            z_cell_index=z_cell_index(end);
            if z_cell_index==N_gridz,
                z_cell_index=z_cell_index-1;
            end
        end
        
        
        
        if j==1,
            cells=[x_cell_index y_cell_index z_cell_index];
        else
            cells=[cells;[x_cell_index y_cell_index z_cell_index]]; %#ok<AGROW>
        end
        
        %display([x_cell_index y_cell_index z_cell_index]);
        MART_coefficient_matrix_3index(i,x_cell_index,y_cell_index,z_cell_index)=Segment_Lengths(j);
    end
    
    [cells_new,~,~] = unique(cells,'rows','stable');
    
    if length(cells_new)~=length(cells),
        fprintf('Ray %d cell 判断错误\n', i);
    else
        fprintf('Ray %d cell 判断正确\n', i);
        clear cells_new;
    end
    
    cells_Rays{i,1}=cells;
    cells=[];
    
    MART_coefficient_matrix(i,:)=reshape(MART_coefficient_matrix_3index(i,:,:,:),1,N_net);
    
end

Grid_raysinfo.Grid_intersection_Rays=Grid_intersection_Rays;
Grid_raysinfo.Segment_Lengths_Rays=Segment_Lengths_Rays;
Grid_raysinfo.cells_Rays=cells_Rays;
Grid_raysinfo.MART_coefficient_matrix_3index=MART_coefficient_matrix_3index;
Grid_raysinfo.Grid_intersection_sp=Grid_intersection_sp;


MART_coefficient_matrix=MART_coefficient_matrix/1.0e3; % m =====> km
% whos  Grid_intersection_Rays Segment_Lengths_Rays cells_Rays MART_coefficient_matrix_3index MART_coefficient_matrix

% save('Ion_MART_coefficient_matrix.mat','MART_coefficient_matrix');



end


