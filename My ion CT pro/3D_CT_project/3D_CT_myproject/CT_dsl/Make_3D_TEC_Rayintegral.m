%% set up initial parameter
clear all
load('MART_coefficient_matrix.mat');

[Nrays,N_net]=size(MART_coefficient_matrix);
% N_net=2080; %电子密度网格数量
% Nrays=1200; %射线数量

A=MART_coefficient_matrix*1000;  %国际单位米
load('Ne2.mat');% 电离层经验模型
load('Ne22.mat');%加上扰动的实际的电离层模型

iri_Ne=repmat(Ne2,10,1);
Ne_initial=reshape(iri_Ne,N_net,1);%初始值 初始估计所有电子密度网格分布
X0=Ne_initial;

Ne_true=repmat(Ne22,10,1);
Ne_true=reshape(Ne_true,N_net,1);  %实际电子密度分布

TEC_R=A*Ne_initial;

save('3D_TEC_R.mat','TEC_R');










