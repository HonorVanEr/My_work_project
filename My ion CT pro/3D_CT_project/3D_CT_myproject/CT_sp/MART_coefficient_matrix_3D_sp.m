function [Grid_raysinfo, MART_coefficient_matrix] = MART_coefficient_matrix_3D_sp(Rays,Xg,Yg,Zg)

%����ֱ������ϵ�¶���������ռ���ά����Ľ����Լ�ÿ������Ľؾ�
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
%   Rays Nx2X3     N ���ߵ�����
%   Rays(i, 1, :)  ��i�����ߵ�������� ̨վλ��
%   Rays(i, 2, :)  ��i�����ߵ��յ����� ����λ��
% �������
%   Xg    lat γ�ȸ��  γ�ȷ�Χ -90����γ�� �� 90����γ��deg
%   Yg    lon ���ȸ��  ���ȷ�Χ ����-180 �� ����180 deg
%   Zg    �߶� �߶ȸ��  ����ĸ߶ȷ�Χ ����ĵľ��� km
%
% OUTPUTS:
%      Grid_raysinfo �ṹ��
%        Grid_i ntersection_Rays ��������ֱ������ϵ
%        Grid_intersection_sp   ��������������ϵ
%        Segment_Lengths_Rays ���߽�ȡ���߶νؾ�
%        cells_Rays    �ؾ����ά�±�
%        MART_coefficient_matrix_3index ϵ������ N*Nx*Ny*Nz



%      MART_coefficient_matrix ϵ������   Nx(Nx*Ny*Nz) ��λkm
%
% Examples:
% Xg=[14.5:0.1:15.5];
% Yg=[70:1:80];
% Hight=[100:20:500];%�������ĸ߶� km
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
Re=6371.2;    %����뾶�����̨վ��λ��


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
    %�߶�ƽ�潻��
    [x_height,y_height,z_height]=Ion_line_plane_intersect(SC_point,GS_point,'height',height);
    xy_intersection=[x_height;y_height;z_height]';
    if size(xy_intersection,2) ~= 3,
        xy_intersection=reshape(xy_intersection,length(xy_intersection)/3,3);
    end
    
    %����
    %     [az,el,r] = cart2sph(xy_intersection(:,1),xy_intersection(:,2),xy_intersection(:,3));
    %     lat_test=rad2deg(el); lon_test=rad2deg(az);
    
    if isempty(xy_intersection),xy_intersection=[];end
    
    %����ƽ�潻��
    [x_lon,y_lon,z_lon]=Ion_line_plane_intersect(SC_point,GS_point,'lon',Yg); %#ok<*ASGLU>
    xz_intersection=[x_lon;y_lon;z_lon]';
    if size(xz_intersection,2) ~= 3,
        xz_intersection=reshape(xz_intersection,length(xz_intersection)/3,3);
    end
    
    if isempty(xz_intersection),xz_intersection=[];end
    
    %γ��ƽ�潻��
    [x_lat,y_lat,z_lat]=Ion_line_plane_intersect(SC_point,GS_point,'lat',Xg);
    yz_intersection=[x_lat;y_lat;z_lat]';
    if size(yz_intersection,2) ~= 3,
        yz_intersection=reshape(yz_intersection,length(yz_intersection)/3,3);
    end
    
    
    if isempty(yz_intersection),yz_intersection=[];end
    
    
    %���н����������������ϰ�һ��ά�����м��� �ų���ͬ����
    Grid_intersection=[xy_intersection;xz_intersection;yz_intersection];
    [Grid_intersection, ~, ~] = unique(Grid_intersection,'rows','stable');    %�Զ����������ˣ�
    
    %�޳�nanֵ
    Grid_intersection(any(isnan(Grid_intersection),2),:)=[];
    
    if isempty(Grid_intersection),
        continue
    end
    
    
    %��X��������
    if line_vet(1)>0,
        temp=sortrows(Grid_intersection,1);
    else
        temp=sortrows(Grid_intersection,-1);
    end
    %��Y��������
    if line_vet(2)>0,
        temp2=sortrows(temp,2);
    else
        temp2=sortrows(temp,-2);
    end
    %��Z��������
    if line_vet(3)>0,
        temp3=sortrows(temp2,3);
    else
        temp3=sortrows(temp2,-3);
    end
    
    %break
    Grid_intersection=temp3;
    clear temp temp2 temp3
    
    
    %�ҵ�λ�������ϵĵ�
    x_ind=find(Grid_intersection(:,1)>=min([x0 x1]) & Grid_intersection(:,1)<=max([x0 x1]));
    y_ind=find(Grid_intersection(:,2)>=min([y0 y1]) & Grid_intersection(:,2)<=max([y0 y1]));
    z_ind=find(Grid_intersection(:,3)>=min([z0 z1]) & Grid_intersection(:,3)<=max([z0 z1]));
    Grid_intersection=Grid_intersection(intersect(intersect(x_ind,y_ind),z_ind),:);
    
    %�ҵ�λ�������ڵĵ� ����ֱ������ϵ����� ��γ������ϵ����
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
    
    
    %����ǳ�������Ϊ��һ����
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
    
    
    %����γ�Ⱦ��Ȳ����ϸ���������
    for j=1:length(Segment_Lengths);
        %�ж��������е�λ������cell ���ַ���
        
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
            
            %���Ȳ�������������
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
        fprintf('Ray %d cell �жϴ���\n', i);
    else
        fprintf('Ray %d cell �ж���ȷ\n', i);
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


