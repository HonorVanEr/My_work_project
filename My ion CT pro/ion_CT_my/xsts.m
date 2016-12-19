clc
clear all

%%satlliate parameters
G=6.67*10^(-12); %万有引力常量
Mz=5.976*10^24;   %地球质量 kg
T=2*pi*(6871200^1.5)/(G*Mz)^0.5;%卫星运动周期开普勒定律
f=1/T*360;               %经度信息
sat=14:f:15.5;           %方位角信息 1x75 个值走过75个网格
Sat=sat*pi/180;          %弧度制
Rs=6871.2;               %卫星离地心的距离


%%receive  parameters81.7459
% rec=[-5  -4 -3  -2 -1 0  1 2 3 4 5];
rec=14.6:0.02:14.9; %布了length(rec)个站
rec=sort(rec);
Rec=rec*pi/180;  %变成弧度制
Rr=6371.2;       %地球半径

%%net paramters
dw=0.0288;               %网格经度间隔
Hight=(Rr+100):10:(Rr+500); %网格的高度范围 1X41 40
Lat=14:dw:15.5;             %网格的方位网格 1X53 52
M=length(rec)*length(Sat); %射线数目 1200
N=(length(Hight)-1)*(length(Lat)-1); %网格数目 40x50=2080
LMN=zeros(M,N);                      %矩阵 1200x2080



jm=1;                                %第一条射线经过的网格点数截距测量
for i=1:length(sat);                 %i=1代表卫星第一个运动方位角位置下标
    for j=1:length(rec);             %j=1代表第一个台站的方位角位置下标
        R=[];
        k=(Rs*sin(Sat(i))-Rr*sin(Rec(j)))/(Rs*cos(Sat(i))-Rr*cos(Rec(j)));%方程斜率
        radial=sqrt(Rs^2+Rr^2-2*Rs*Rr*cos(Sat(i)-Rec(j))); %卫星与台站连线距离
        alfa=acos((Rs^2+radial^2-Rr^2)/(2*Rs*radial));     %卫星与台站连线和卫星地心连线夹角
        %卫星在接收台站的右侧 找到两者之间方位方向的网格点位置
        %找到卫星离得最近左侧的网格
        %找到台站离得最近的右侧的网格
        if sat(i)>=rec(j)
            sita1=ceil(rec(j)/dw);
            sita1=sita1*dw;
            sita2=floor(sat(i)/dw);
            sita2=sita2*dw;
            sita=sita1:dw:sita2;
            sita=sort(sita,'descend');
            Sita=sita*pi/180;
            Sita=sort(Sita,'descend');
            
            %计算地球与经过这些网格射线上点的径向距离
            for q=1:length(sita)
                R(q)=(k*Rr*cos(Rec(j))-Rr*sin(Rec(j)))/(k*cos(Sita(q))-sin(Sita(q)));%radial square
                R=R(R>Hight(1));
                
            end
            
            %不同高度上的网格对应得方位角 
            for r=1:length(Hight)
                gama(r)=asin(sin(alfa)*Rs/Hight(r));
                gama=gama*180/pi;
                sitaz(r)=alfa*180/pi-gama(r)+sat(i);
            end
            %卫星在接收台站的左侧
            %找到卫星离得最近右侧的网格
            %找到台站离得最近的左侧的网格
        else
            sita1=floor(rec(j)/dw);
            sita1=sita1*dw;
            sita2=ceil(sat(i)/dw);
            sita2=sita2*dw;
            sita=sita2:dw:sita1;
            sita=sort(sita,'descend');
            Sita=sita*pi/180;
            Sita=sort(Sita,'descend');
            
            %计算地球与经过这些网格射线上点的径向距离
            for q=1:length(sita)
                R(q)=(k*Rr*cos(Rec(j))-Rr*sin(Rec(j)))/(k*cos(Sita(q))-sin(Sita(q)));%radial square
            end
            R=R(R>Hight(1));%超过Hight(1)才是有效的
            
            
            %不同高度上的网格对应得方位角 不太确定？
            for r=1:length(Hight)
                gama(r)=asin(sin(alfa)*Rs/Hight(r));
                gama(r)=gama(r)*180/pi;
                sitaz(r)=gama(r)-alfa*180/pi+sat(i);
            end
        end
        %卫星和台站方位方向得到的网格点和高度方向得到的点的交集
        sitau=union(sita,sitaz);
        sitau=sort(sitau);
        
        if sat(i)>=rec(j)
            sitau=sitau(sitau>=sitaz(1));
        else
            sitau=sitau(sitau<=sitaz(1));
        end
        tmp=diff(sitau);
        for l=i:length(tmp)
            if tmp(l)<1e-10
                sitau(l)=[];
            end
        end
        if sat(i)>=rec(j)
            sitau=sort(sitau) ;
        else
            sitau=sort(sitau,'descend') ;
        end
        Sitau=sitau*pi/180;
        Ru=union(R,Hight);
        Ru=sort(Ru);
        %该射线穿过的有效方位角和高度方向
        for t=1:length(Sitau)-1
            L=(Ru(t+1)*sin(Sitau(t+1))-Ru(t)*sin(Sitau(t)))^2+(Ru(t+1)*cos(Sitau(t+1))-Ru(t)*cos(Sitau(t)))^2; %求每个网格的截距
            L=sqrt(L);
            m=floor((Ru(t)-Hight(1))/10)+1;       %高度上第几个网格点下标
            tmp=(sitau(t)-Lat(1))/dw;             %纬度方向第几个网点的下标
            if mod(tmp,1)==0
                if sat(i)>=rec(j)                 %卫星在台站右侧下标增加一个
                    n=floor((sitau(t)-Lat(1))/dw)+1;
                else
                    n=floor((sitau(t)-Lat(1))/dw); %卫星在台站左侧下标不变
                end
            else
                n=floor((sitau(t)-Lat(1))/dw)+1;
            end
            
            jn=((n-1)*(length(Hight)-1)+m);    %找到穿过的2080个网格中网格点的一维下标
            LMN(jm,jn)=L;%系数矩阵
            
        end
        %        软件测试：
        %        tu=LMN(jm,:);
        %        tu=reshape(tu,40,52);
        %        spy(tu)
        jm=jm+1;
    end
end
tu=LMN(80,:);
tu=reshape(tu,40,52);
spy(tu)
save('xishu.mat','LMN');

