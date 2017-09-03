function [Grid_raysinfo, MART_coefficient_matrix] = MART_coefficient_matrix_3D(Rays,Xg,Yg,Zg)

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
%   Xg
%   Yg
%   Zg
%
% OUTPUTS:
%      Grid_raysinfo �ṹ��
%      MART_coefficient_matrix ϵ������   Nx(Nx*Ny*Nz)
%
% Examples:
%
% Rays(1,1,:)=[-5;1;21];
% Rays(1,2,:)=[122;79;69];
%
% Xg = 10*[-2 1 3 5 6 8 10 13];
% Yg = 10*[-1 1 3 5 6 8 10];
% Zg = 10*[ 0 1 3 5 6 8];
% [Grid_raysinfo, MART_coefficient_matrix]=MART_coefficient_matrix_3D(Rays,Xg,Yg,Zg);
%





% MART coefficient_matrix
%�±��Ȱ�X�����ٰ�Y�������Z����
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
    
    %xy ƽ�潻��
    planexy_n=[0 0 1];planexy_n=repmat(planexy_n,N_gridz,1);
    plane_point=zeros(N_gridz,3); plane_point(:,3)=Zg;
    [x,y,z]=find_intersection(planexy_n,plane_point,line_vet,line_point);
    xy_intersection=[x,y,z];
    
    %xz ƽ�潻��
    planexz_n=[0 1 0];planexz_n=repmat(planexz_n,N_gridy,1);
    plane_point=zeros(N_gridy,3); plane_point(:,2)=Yg;
    [x,y,z]=find_intersection(planexz_n,plane_point,line_vet,line_point);
    xz_intersection=[x,y,z];
    
    %yz ƽ�潻��
    planeyz_n=[1 0 0];planeyz_n=repmat(planeyz_n,N_gridx,1);
    plane_point=zeros(N_gridx,3); plane_point(:,1)=Xg;
    [x,y,z]=find_intersection(planeyz_n,plane_point,line_vet,line_point);
    yz_intersection=[x,y,z];
    
    %���н����������������ϰ�һ��ά�����м��� �ų���ͬ����
    Grid_intersection=[xy_intersection;xz_intersection;yz_intersection];
    [Grid_intersection, ~, ~] = unique(Grid_intersection,'rows','stable');
    
    %�޳�nanֵ
    Grid_intersection(any(isnan(Grid_intersection),2),:)=[];
    
    if isempty(Grid_intersection),
        continue
    end
    
    
    %����ǳ�������Ϊ��һ����
    tol=eps;
    
    [Grid_intersection,~,~] = uniquetol(Grid_intersection,tol,'ByRows',true);
    
    
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
    
    %�ҵ�λ�������ڵĵ�
    x_ind=find(Grid_intersection(:,1)>=Xgrid_Min & Grid_intersection(:,1)<=Xgrid_Max);
    y_ind=find(Grid_intersection(:,2)>=Ygrid_Min & Grid_intersection(:,2)<=Ygrid_Max);
    z_ind=find(Grid_intersection(:,3)>=Zgrid_Min & Grid_intersection(:,3)<=Zgrid_Max);
    
    Grid_intersection=Grid_intersection(intersect(intersect(x_ind,y_ind),z_ind),:);
    
    %     %�ж�����ǳ��ӽ��ĵ�
    %     tmp=diff(Grid_intersection, 1);
    %     for l=1:length(tmp)
    %         if tmp(l)<1e-10
    %             Grid_intersection(l)=[];
    %         end
    %     end
    
    
    Grid_intersection_Rays{i,1}=Grid_intersection;
    
    
    
    if (size(Grid_intersection,1)>1),
        dRi = diff(Grid_intersection, 1);
        Segment_Lengths = sqrt(sum(dRi.^2, 2));
        Segment_Lengths_Rays{i,1}=Segment_Lengths; %#ok<SAGROW>
    else
        continue
    end
    
    for j=1:length(Segment_Lengths);
        %�ж��������е�λ������cell
        x_cell_index=find(Xg<=(Grid_intersection(j,1)+Grid_intersection(j+1,1))/2);
        x_cell_index=x_cell_index(end);
        if x_cell_index==N_gridx,
            x_cell_index=x_cell_index-1;
        end
        
        y_cell_index=find(Yg<=(Grid_intersection(j,2)+Grid_intersection(j+1,2))/2);
        y_cell_index=y_cell_index(end);
        if y_cell_index==N_gridy,
            y_cell_index=y_cell_index-1;
        end
        
        z_cell_index=find(Zg<=(Grid_intersection(j,3)+Grid_intersection(j+1,3))/2);
        z_cell_index=z_cell_index(end);
        if z_cell_index==N_gridz,
            z_cell_index=z_cell_index-1;
        end
        
        if j==1,
            cells=[x_cell_index y_cell_index z_cell_index];
        else
            cells=[cells;[x_cell_index y_cell_index z_cell_index]]; %#ok<AGROW>
        end
        
        %display([x_cell_index y_cell_index z_cell_index]);
        MART_coefficient_matrix_3index(i,x_cell_index,y_cell_index,z_cell_index)=Segment_Lengths(j);
        
    end
    cells_Rays{i,1}=cells;
    [C, ~, ~] = unique(cells,'rows');
    
    if size(C,1)~=size(cells,1),
        fprintf('Ray %d cell �ظ�\n', i);
    end
    
    MART_coefficient_matrix(i,:)=reshape(MART_coefficient_matrix_3index(i,:,:,:),1,N_net);
    
end
Grid_raysinfo.Grid_intersection_Rays=Grid_intersection_Rays;
Grid_raysinfo.Segment_Lengths_Rays=Segment_Lengths_Rays;
Grid_raysinfo.cells_Rays=cells_Rays;
Grid_raysinfo.MART_coefficient_matrix_3index=MART_coefficient_matrix_3index;

%whos  Grid_intersection_Rays Segment_Lengths_Rays cells_Rays MART_coefficient_matrix_3index MART_coefficient_matrix

end


