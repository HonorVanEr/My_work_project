%% set up initial parameter
clc
clear all
close all
L=2080; %�����ܶ���������
I=1200; %��������
lmt=1;  %�ɳ�����



load('xishu.mat');
B=LMN;
B=B*1000;  %����
load('Ne2.mat');% ����㾭��ģ��
load('Ne22.mat');%�����Ŷ���ʵ�ʵĵ����ģ��

TEC_initial=reshape(Ne2,L,1);%��ʼֵ ��ʼ�������е����ܶ�����ֲ�
X2=TEC_initial;

%% begin run
tic

x=reshape(Ne22,L,1);  %ʵ�ʵ����ܶȷֲ�


load TEC_R.txt %����SAR���ݵõ���TECֵ
yyy=TEC_R;

yyy=reshape(yyy,75,16);% ���ǵ����� 16��̨վ 75��ʾ���ǹ�� Azimuth Direcction�߹�������
yyy=yyy.';
yyy=reshape(yyy,I,1);


for j=1:3  %û���10�ε����������3�ε���
    for i=1:I
        cha2=yyy(i)/dot(X2,B(i,:));  %y(i)��Y(i)�ֱ���������Ƶ����˫Ƶ����������TECʵ��ֵ������Y(i)�Ǵ�������
        mo=norm(B(i,:));
        for l=1:L
            X2(l)=X2(l)*(cha2^(lmt*B(i,l)'/mo));
        end
    end
    j
end

X2=reshape(X2,40,52);   %���ݵõ��ĵ����ܶȷֲ�
x=reshape(x,40,52);     %%�����Ŷ���ʵ�ʵĵ����ģ��

v=0:.12e11:1.56e11;
Azi=0:0.0288*82.4891:1.5*82.4891-0.0288*82.4891;
figure

contourf(Azi(9:44)-Azi(9),205:5:400,x(:,9:44),v);
h=colorbar;
set(get(h,'title'),'string');
axis([0 80 205 400])
xlabel('Azimuth Direcction/km');
ylabel('Altitude/km');

figure
contourf(Azi(9:44)-Azi(9),205:5:400,X2(:,9:44),v);
h=colorbar;
set(get(h,'title'),'string');
axis([0 80 205 400])
xlabel('Azimuth Direcction/km');
ylabel('Altitude/km');

save('x.mat','x');
save('X2.mat','X2')
toc