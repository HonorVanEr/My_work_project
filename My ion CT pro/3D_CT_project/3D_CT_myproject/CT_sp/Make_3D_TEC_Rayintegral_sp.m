%% set up initial parameter
clear all
load('MART_coefficient_matrix_sp.mat');

[Nrays,N_net]=size(MART_coefficient_matrix);
% N_net=2080; %�����ܶ���������
% Nrays=1200; %��������

A=MART_coefficient_matrix*1000;  %���ʵ�λ��
load('Ne2.mat');% ����㾭��ģ��
load('Ne22.mat');%�����Ŷ���ʵ�ʵĵ����ģ��

Ne2=Ne2';
iri_Ne=repmat(Ne2,10,1);
iri_Ne_3d=reshape(iri_Ne,Nx,Ny,Nz);  %����3ά����ʽ
Ne_initial=reshape(iri_Ne_3d,N_net,1);%��ʼֵ ��ʼ�������е����ܶ�����ֲ�
X0=Ne_initial;


Ne22=Ne22';
Ne_true=repmat(Ne22,10,1);
Ne_true_3d=reshape(Ne_true,Nx,Ny,Nz);
Ne_true=reshape(Ne_true,N_net,1);  %ʵ�ʵ����ܶȷֲ�

save('3D_TEC_R_sp.mat','TEC_R');










