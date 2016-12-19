clc
clear all

%%satlliate parameters
G=6.67*10^(-12); %������������
Mz=5.976*10^24;   %�������� kg
T=2*pi*(6871200^1.5)/(G*Mz)^0.5;%�����˶����ڿ����ն���
f=1/T*360;               %������Ϣ
sat=14:f:15.5;           %��λ����Ϣ 1x75 ��ֵ�߹�75������
Sat=sat*pi/180;          %������
Rs=6871.2;               %��������ĵľ���


%%receive  parameters81.7459
% rec=[-5  -4 -3  -2 -1 0  1 2 3 4 5];
rec=14.6:0.02:14.9; %����length(rec)��վ
rec=sort(rec);
Rec=rec*pi/180;  %��ɻ�����
Rr=6371.2;       %����뾶

%%net paramters
dw=0.0288;               %���񾭶ȼ��
Hight=(Rr+100):10:(Rr+500); %����ĸ߶ȷ�Χ 1X41 40
Lat=14:dw:15.5;             %����ķ�λ���� 1X53 52
M=length(rec)*length(Sat); %������Ŀ 1200
N=(length(Hight)-1)*(length(Lat)-1); %������Ŀ 40x50=2080
LMN=zeros(M,N);                      %���� 1200x2080



jm=1;                                %��һ�����߾�������������ؾ����
for i=1:length(sat);                 %i=1�������ǵ�һ���˶���λ��λ���±�
    for j=1:length(rec);             %j=1������һ��̨վ�ķ�λ��λ���±�
        R=[];
        k=(Rs*sin(Sat(i))-Rr*sin(Rec(j)))/(Rs*cos(Sat(i))-Rr*cos(Rec(j)));%����б��
        radial=sqrt(Rs^2+Rr^2-2*Rs*Rr*cos(Sat(i)-Rec(j))); %������̨վ���߾���
        alfa=acos((Rs^2+radial^2-Rr^2)/(2*Rs*radial));     %������̨վ���ߺ����ǵ������߼н�
        %�����ڽ���̨վ���Ҳ� �ҵ�����֮�䷽λ����������λ��
        %�ҵ�������������������
        %�ҵ�̨վ���������Ҳ������
        if sat(i)>=rec(j)
            sita1=ceil(rec(j)/dw);
            sita1=sita1*dw;
            sita2=floor(sat(i)/dw);
            sita2=sita2*dw;
            sita=sita1:dw:sita2;
            sita=sort(sita,'descend');
            Sita=sita*pi/180;
            Sita=sort(Sita,'descend');
            
            %��������뾭����Щ���������ϵ�ľ������
            for q=1:length(sita)
                R(q)=(k*Rr*cos(Rec(j))-Rr*sin(Rec(j)))/(k*cos(Sita(q))-sin(Sita(q)));%radial square
                R=R(R>Hight(1));
                
            end
            
            %��ͬ�߶��ϵ������Ӧ�÷�λ�� 
            for r=1:length(Hight)
                gama(r)=asin(sin(alfa)*Rs/Hight(r));
                gama=gama*180/pi;
                sitaz(r)=alfa*180/pi-gama(r)+sat(i);
            end
            %�����ڽ���̨վ�����
            %�ҵ������������Ҳ������
            %�ҵ�̨վ����������������
        else
            sita1=floor(rec(j)/dw);
            sita1=sita1*dw;
            sita2=ceil(sat(i)/dw);
            sita2=sita2*dw;
            sita=sita2:dw:sita1;
            sita=sort(sita,'descend');
            Sita=sita*pi/180;
            Sita=sort(Sita,'descend');
            
            %��������뾭����Щ���������ϵ�ľ������
            for q=1:length(sita)
                R(q)=(k*Rr*cos(Rec(j))-Rr*sin(Rec(j)))/(k*cos(Sita(q))-sin(Sita(q)));%radial square
            end
            R=R(R>Hight(1));%����Hight(1)������Ч��
            
            
            %��ͬ�߶��ϵ������Ӧ�÷�λ�� ��̫ȷ����
            for r=1:length(Hight)
                gama(r)=asin(sin(alfa)*Rs/Hight(r));
                gama(r)=gama(r)*180/pi;
                sitaz(r)=gama(r)-alfa*180/pi+sat(i);
            end
        end
        %���Ǻ�̨վ��λ����õ��������͸߶ȷ���õ��ĵ�Ľ���
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
        %�����ߴ�������Ч��λ�Ǻ͸߶ȷ���
        for t=1:length(Sitau)-1
            L=(Ru(t+1)*sin(Sitau(t+1))-Ru(t)*sin(Sitau(t)))^2+(Ru(t+1)*cos(Sitau(t+1))-Ru(t)*cos(Sitau(t)))^2; %��ÿ������Ľؾ�
            L=sqrt(L);
            m=floor((Ru(t)-Hight(1))/10)+1;       %�߶��ϵڼ���������±�
            tmp=(sitau(t)-Lat(1))/dw;             %γ�ȷ���ڼ���������±�
            if mod(tmp,1)==0
                if sat(i)>=rec(j)                 %������̨վ�Ҳ��±�����һ��
                    n=floor((sitau(t)-Lat(1))/dw)+1;
                else
                    n=floor((sitau(t)-Lat(1))/dw); %������̨վ����±겻��
                end
            else
                n=floor((sitau(t)-Lat(1))/dw)+1;
            end
            
            jn=((n-1)*(length(Hight)-1)+m);    %�ҵ�������2080��������������һά�±�
            LMN(jm,jn)=L;%ϵ������
            
        end
        %        �������ԣ�
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
