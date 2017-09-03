 %% line_plane_intersect
function [x,y,z]=Ion_line_plane_intersect(SC_point,GS_point,varargin)

% Ion_line_plane_intersect - 与高度面 经度面 纬度面的交点
%
%   [x,y,z]=Ion_line_plane_intersect(SC_point,GS_point,varargin)
%     Input:
%       SC_point - 卫星位置1x3 WGS坐标系
%       GS_point - 地面台站位置1x3 WGS坐标系
%       varargin - 'height'   height 高度面 离地球表面的位置 单位 m
%                - 'lon'      lon 与起始子午面的夹角经度 角度制
%                - 'lat'      lat 纬度W 角度制
%                - 'vis'      窗口显示出交点信息
%
%     Output:
%       [x,y,z] WGS坐标系
%
%
%   Examples:
%     [x_height,y_height,z_height]=Ion_line_plane_intersect(SC_point,GS_point,'height',height);
%     [x_lon,y_lon,z_lon]=Ion_line_plane_intersect(SC_point,GS_point,'lon',Yg);
%     [x_lat,y_lat,z_lat]=Ion_line_plane_intersect(SC_point,GS_point,'lat',Xg);
%
%[x_height,y_height,z_height]=Ion_line_plane_intersect(SC_point,GS_point,'height',height,'vis');
%end

%initail
SC_x=SC_point(1);
SC_y=SC_point(2);
SC_z=SC_point(3);
GS_x=GS_point(1);
GS_y=GS_point(2);
GS_z=GS_point(3);
Re=6371.2*1000.0; %m

%Re=6377000; %m  地球平均半径

if ~isempty(varargin),
    switch lower(varargin{1})
        case 'height'               % 高度球面
            height=varargin{2};
            x_temp=[];
            y_temp=[];
            z_temp=[];
            
            for i=1:length(height)
                A2=((GS_y - SC_y)^2/(GS_x - SC_x)^2 + (GS_z - SC_z)^2/(GS_x - SC_x)^2 + 1);
                A1=((2*(GS_y - SC_y)*(GS_y - (GS_x*(GS_y - SC_y))/(GS_x - SC_x)))/(GS_x - SC_x) + (2*(GS_z - SC_z)*(GS_z - (GS_x*(GS_z - SC_z))/(GS_x - SC_x)))/(GS_x - SC_x));
                A0=(GS_y - (GS_x*(GS_y - SC_y))/(GS_x - SC_x))^2 - (Re + height(i))^2 + (GS_z - (GS_x*(GS_z - SC_z))/(GS_x - SC_x))^2;
                Croots=roots([A2 A1 A0]);
                Croots=Croots(abs(imag(Croots))<=eps);   %判断是否有复数根并剔除
                
                x=Croots(Croots>=min([SC_x,GS_x]) & Croots<=max([SC_x,GS_x]));
                
                
                y=GS_y+(x-GS_x)*(SC_y-GS_y)/(SC_x-GS_x);
                z=GS_z+(x-GS_x)*(SC_z-GS_z)/(SC_x-GS_x);
                
                x_temp=[x_temp x]; %#ok<AGROW>
                y_temp=[y_temp y]; %#ok<*NOPRT>
                z_temp=[z_temp z]; %#ok<*AGROW>
                
                x_height=x;
                y_height=y;
                z_height=z;
                if length(varargin)>2 && sum(lower(varargin{3})=='vis')/3
                    fprintf('Intersect point Height m: X = %2f, Y  = %2f, Z = %2f\n',x_height,y_height,z_height);
                end
            end
            x=x_temp;
            y=y_temp;
            z=z_temp;
            
        case 'lon'          % 经度面 经度面是一个穿过Z轴的平面
            lon=varargin{2};
            x_temp=[];
            y_temp=[];
            z_temp=[];
            
            for i=1:length(lon)
                if SC_x~=GS_x
                    if (abs(lon(i))==90)
                        x=0;
                        y=GS_y+(x-GS_x)*(SC_y-GS_y)/(SC_x-GS_x);
                        z=GS_z+(x-GS_x)*(SC_z-GS_z)/(SC_x-GS_x);
                    else
                        B1=(tand(lon(i)) - (GS_y - SC_y)/(GS_x - SC_x));
                        B0=(GS_x*(GS_y - SC_y))/(GS_x - SC_x) - GS_y;
                        Croots=-B0/B1;
                        Croots=Croots(abs(imag(Croots))<=eps);   %判断是否有复数根并剔除
                        
                        x=Croots(Croots>=min([SC_x,GS_x]) & Croots<=max([SC_x,GS_x]));
                        y=GS_y+(x-GS_x)*(SC_y-GS_y)/(SC_x-GS_x);
                        z=GS_z+(x-GS_x)*(SC_z-GS_z)/(SC_x-GS_x);
                    end
                else
                    x=[];
                    y=[];
                    z=[];
                end
                x_temp=[x_temp x]; %#ok<AGROW>
                y_temp=[y_temp y]; %#ok<*NOPRT>
                z_temp=[z_temp z]; %#ok<*AGROW>
                
                x_lon=x;
                y_lon=y;
                z_lon=z;
                if length(varargin)>2 && sum(lower(varargin{3})=='vis')/3
                    fprintf('Intersect point LONGITUDE: X = %2f, Y  = %2f, Z = %2f\n',x_lon,y_lon,z_lon);
                end
            end
            x=x_temp;
            y=y_temp;
            z=z_temp;
            
        case 'lat'   % 纬度面 纬度面是由某一点与地心的连线绕Z轴旋转形成的圆锥面
            lat=varargin{2};
            x_temp=[];
            y_temp=[];
            z_temp=[];
            
            for i=1:length(lat)
                
                C2=((GS_z - SC_z)^2/(GS_x - SC_x)^2 - tand(lat(i))^2*((GS_y - SC_y)^2/(GS_x - SC_x)^2 + 1));
                C1=((2*(GS_z - SC_z)*(GS_z - (GS_x*(GS_z - SC_z))/(GS_x - SC_x)))/(GS_x - SC_x) - (2*tand(lat(i))^2*(GS_y - SC_y)*(GS_y - (GS_x*(GS_y - SC_y))/(GS_x - SC_x)))/(GS_x - SC_x));
                C0=(GS_z - (GS_x*(GS_z - SC_z))/(GS_x - SC_x))^2 - tand(lat(i))^2*(GS_y - (GS_x*(GS_y - SC_y))/(GS_x - SC_x))^2;
                Croots=roots([C2 C1 C0]);
                Croots=Croots(abs(imag(Croots))<=eps);   %判断是否有复数根并剔除
                
                x=Croots(Croots>=min([SC_x,GS_x]) & Croots<=max([SC_x,GS_x]));
                y=GS_y+(x-GS_x)*(SC_y-GS_y)/(SC_x-GS_x);
                z=GS_z+(x-GS_x)*(SC_z-GS_z)/(SC_x-GS_x);
                
                x_temp=[x_temp x]; %#ok<AGROW>
                y_temp=[y_temp y]; %#ok<*NOPRT>
                z_temp=[z_temp z]; %#ok<*AGROW>
                
                
                x_lat=x;
                y_lat=y;
                z_lat=z;
                if length(varargin)>2 && sum(lower(varargin{3})=='vis')/3
                    fprintf('Intersect point Latitude: X = %2f, Y  = %2f, Z = %2f\n',x_lat,y_lat,z_lat);
                end
            end
            
            x=x_temp;
            y=y_temp;
            z=z_temp;
            
    end
    
end


end




