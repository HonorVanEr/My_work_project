%% set up initial parameter
clc
clear all
close all
L=2080; %电子密度网格数量
I=1200; %射线数量
lmt=1;  %松弛因子



load('xishu.mat');
B=LMN;
B=B*1000;  %公里
load('Ne2.mat');% 电离层经验模型
load('Ne22.mat');%加上扰动的实际的电离层模型

TEC_initial=reshape(Ne2,L,1);%初始值 初始估计所有电子密度网格分布
X2=TEC_initial;

%% begin run
tic

x=reshape(Ne22,L,1);  %实际电子密度分布


load TEC_R.txt %极化SAR反演得到的TEC值
yyy=TEC_R;

yyy=reshape(yyy,75,16);% 卫星的网格 16个台站 75表示卫星轨道 Azimuth Direcction走过的网格
yyy=yyy.';
yyy=reshape(yyy,I,1);


for j=1:3  %没误差10次迭代，有误差3次迭代
    for i=1:I
        cha2=yyy(i)/dot(X2,B(i,:));  %y(i)与Y(i)分别是利用三频法与双频法测量到的TEC实际值，其中Y(i)是带有误差的
        mo=norm(B(i,:));
        for l=1:L
            X2(l)=X2(l)*(cha2^(lmt*B(i,l)'/mo));
        end
    end
    j
end

X2=reshape(X2,40,52);   %反演得到的电子密度分布
x=reshape(x,40,52);     %%加上扰动的实际的电离层模型

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