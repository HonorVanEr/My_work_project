%% line_plane_intersect
function [x,y,z]=line_plane_intersect(SC_point,GS_point,lat,lon,height)

SC_x=SC_point(1);
SC_y=SC_point(2);
SC_z=SC_point(3);
GS_x=GS_point(1);
GS_y=GS_point(2);
GS_z=GS_point(3);
Re=6371;

% SC_x=1;
% SC_y=1;
% SC_z=2;
% GS_x=3;
% GS_y=5;
% GS_z=8;
% Re=1;
% height=2;
% lat=50;
% lon=80;

% 高度球面
A2=((GS_y - SC_y)^2/(GS_x - SC_x)^2 + (GS_z - SC_z)^2/(GS_x - SC_x)^2 + 1);
A1=((2*(GS_y - SC_y)*(GS_y - (GS_x*(GS_y - SC_y))/(GS_x - SC_x)))/(GS_x - SC_x) + (2*(GS_z - SC_z)*(GS_z - (GS_x*(GS_z - SC_z))/(GS_x - SC_x)))/(GS_x - SC_x));
A0=(GS_y - (GS_x*(GS_y - SC_y))/(GS_x - SC_x))^2 - (Re + height)^2 + (GS_z - (GS_x*(GS_z - SC_z))/(GS_x - SC_x))^2;
Croots=roots([A2 A1 A0]);
x=Croots(Croots>=min([SC_x,GS_x]) & Croots<=max([SC_x,GS_x]));
y=GS_y+(x-GS_x)*(SC_y-GS_y)/(SC_x-GS_x);
z=GS_z+(x-GS_x)*(SC_z-GS_z)/(SC_x-GS_x);
x_height=x;
y_height=y;
z_height=z;

fprintf('Intersect point Height: X = %2f, Y  = %2f, Z = %2f\n',x_height,y_height,z_height);


% 经度面
if SC_x~=GS_x
    if (abs(lon)==90)
        x=0;
        y=GS_y+(x-GS_x)*(SC_y-GS_y)/(SC_x-GS_x);
        z=GS_z+(x-GS_x)*(SC_z-GS_z)/(SC_x-GS_x);
    else
        B1=(tand(lon) - (GS_y - SC_y)/(GS_x - SC_x));
        B0=(GS_x*(GS_y - SC_y))/(GS_x - SC_x) - GS_y;
        Croots=-B0/B1;
        x=Croots(Croots>=min([SC_x,GS_x]) & Croots<=max([SC_x,GS_x]));
        y=GS_y+(x-GS_x)*(SC_y-GS_y)/(SC_x-GS_x);
        z=GS_z+(x-GS_x)*(SC_z-GS_z)/(SC_x-GS_x);
    end
else
    x=[];
    y=[];
    z=[];
end

x_lon=x;
y_lon=y;
z_lon=z;

fprintf('Intersect point LONGITUDE: X = %2f, Y  = %2f, Z = %2f\n',x_lon,y_lon,z_lon);


% 纬度面
C2=((GS_z - SC_z)^2/(GS_x - SC_x)^2 - tand(lat)^2*((GS_y - SC_y)^2/(GS_x - SC_x)^2 + 1));
C1=((2*(GS_z - SC_z)*(GS_z - (GS_x*(GS_z - SC_z))/(GS_x - SC_x)))/(GS_x - SC_x) - (2*tand(lat)^2*(GS_y - SC_y)*(GS_y - (GS_x*(GS_y - SC_y))/(GS_x - SC_x)))/(GS_x - SC_x));
C0=(GS_z - (GS_x*(GS_z - SC_z))/(GS_x - SC_x))^2 - tand(lat)^2*(GS_y - (GS_x*(GS_y - SC_y))/(GS_x - SC_x))^2;
Croots=roots([C2 C1 C0]);
x=Croots(Croots>=min([SC_x,GS_x]) & Croots<=max([SC_x,GS_x]));
y=GS_y+(x-GS_x)*(SC_y-GS_y)/(SC_x-GS_x);
z=GS_z+(x-GS_x)*(SC_z-GS_z)/(SC_x-GS_x);

x_lat=x;
y_lat=y;
z_lat=z;

fprintf('Intersect point Latitude: X = %2f, Y  = %2f, Z = %2f\n',x_lat,y_lat,z_lat);


end


% syms x SC_x SC_y SC_z GS_x GS_y GS_z lat lon height Re
% y=GS_y+(x-GS_x)*(SC_y-GS_y)/(SC_x-GS_x);
% z=GS_z+(x-GS_x)*(SC_z-GS_z)/(SC_x-GS_x);
% EQ1=x^2+y^2+z^2-(Re+height)^2;
% EQ2=tand(lon)*x-y;
% EQ3=z^2-tand(lat)^2*(x^2+y^2);
% collect(EQ1)

