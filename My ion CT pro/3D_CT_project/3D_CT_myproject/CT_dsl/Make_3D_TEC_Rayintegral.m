%% set up initial parameter
clear all
load('MART_coefficient_matrix.mat');

[Nrays,N_net]=size(MART_coefficient_matrix);
% N_net=2080; %�����ܶ���������
% Nrays=1200; %��������

A=MART_coefficient_matrix*1000;  %���ʵ�λ��
load('Ne2.mat');% ����㾭��ģ��
load('Ne22.mat');%�����Ŷ���ʵ�ʵĵ����ģ��

iri_Ne=repmat(Ne2,10,1);
Ne_initial=reshape(iri_Ne,N_net,1);%��ʼֵ ��ʼ�������е����ܶ�����ֲ�
X0=Ne_initial;

Ne_true=repmat(Ne22,10,1);
Ne_true=reshape(Ne_true,N_net,1);  %ʵ�ʵ����ܶȷֲ�

TEC_R=A*Ne_initial;

save('3D_TEC_R.mat','TEC_R');










