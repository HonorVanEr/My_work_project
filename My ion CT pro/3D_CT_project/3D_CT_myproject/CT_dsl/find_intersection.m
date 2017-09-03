
function [x,y,z]=find_intersection(plane_n,plane_point,line_vet,line_point)
%% find_intersection
% Ion_line_plane_intersect - 直角坐标系下与XY XZ YZ 平面交点
%
%   [x,y,z]=find_intersection(plane_n,plane_point,line_vet,line_point)
%     Input:
%       plane_n - 平面法向量 1x3 或者 Nx3
%       plane_point - 平面内的一点 1x3 或者 Nx3
%       line_vet      直线的方向向量
%       line_point    直线通过的点
%
%     Output:
%    [x,y,z] 交点的三维下标
%
%   Examples:
%       plane_n=[0 0 1];  
%       plane_point=[0 0 1];
%       line_vet=[1 1 1];
%       line_point=[0 0 0];%一个面
%       [x,y,z]=find_intersection(plane_n,plane_point,line_vet,line_point)



%多个面
% plane_n=[0 0 1; 0 0 1];
% plane_point=[0 0 1; 0 0 2];
% line_vet=[1 1 1; 1 1 1];
% line_point=[0 0 0; 0 0 1];
% [x,y,z]=find_intersection(plane_n,plane_point,line_vet,line_point)
%



plane_dims=size(plane_n);
line_dims=size(line_vet);
line_pt_dims=size(line_point);
if line_dims(1)~=plane_dims(1)
    line_vet=repmat(line_vet,plane_dims(1),1);
end
if line_pt_dims(1)~=plane_dims(1)
    line_point=repmat(line_point,plane_dims(1),1);
end
ml=line_vet(:,1);
nl=line_vet(:,2);
pl=line_vet(:,3);
x0=line_point(:,1);
y0=line_point(:,2);
z0=line_point(:,3);
x1=plane_point(:,1);
y1=plane_point(:,2);
z1=plane_point(:,3);

ms=plane_n(:,1);
ns=plane_n(:,2);
ps=plane_n(:,3);

B=((ms + (nl.*ns)./ml + (pl.*ps)./ml));
x=-(- ms.*x1 - ns.*(y1 - y0 + (nl.*x0)./ml) - ps.*(z1 - z0 + (pl.*x0)./ml))./B;
y=y0+(nl.*(x-x0))./ml;
z=z0+(pl.*(x-x0))./ml;


%% 解方程步骤
% npts=length(plane_n(:,1));
% x=nan(npts,length(x1));
% y=nan(npts,length(x1));
% z=nan(npts,length(x1));
% for ii=1:npts
%     ms=plane_n(ii,1);
%     ns=plane_n(ii,2);
%     ps=plane_n(ii,3);
%
%     x(ii,:)=-(- ms*x1 - ns*(y1 - y0 + (nl*x0)/ml) - ps*(z1 - z0 + (pl*x0)/ml))/((ms + (nl*ns)/ml + (pl*ps)/ml));
%     y(ii,:)=y0+nl*(x(ii,:)-x0)/ml;
%     z(ii,:)=z0+pl*(x(ii,:)-x0)/ml;
% end
% ml=1;
% nl=1;
% pl=1;
% x0=0;
% y0=0;
% z0=0;
% ms=0;
% ns=0;
% ps=1;
% x1=[0 0 0];
% y1=[0 0 0];
% z1=[1 1 1];
% syms ml nl pl ms ns ps x0 y0 z0 x1 y1 z1 x
%
% y=y0+nl*(x-x0)/ml;
% z=z0+pl*(x-x0)/ml;
% eq4=ms*(x-x1)+ns*(y-y1)+ps*(z-z1);

end


