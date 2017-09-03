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

if 0,
    load TEC_R.txt %极化SAR反演得到的TEC值
    Y=TEC_R;
    Y=reshape(Y,75,16);%  射线数
    Y=Y.';
    Y=reshape(Y,Nrays,1);
else
   load('3D_TEC_R.mat');
   Y=TEC_R;
end


if length(Y)==Nrays,
    display('Y与射线数一样')
end


%mart(p, A, x0, relax=1., iters=1, tol=0):
%   """Performs multiplicitive algebraic reconstruction (MART) of image`x` given
%     projection data `Y` and a projection matrix `A` which satisfies
%     :math:`\vec{Y} = A\vec{s}`

%     Notes: must assume positivity of image--which makes sense for ionosphere tomography, thus x0 must be > 0
%     x_j^{k+1} = x_j^k * (y_i / (\sum_j a_{ij}x_j^k))^{\gamma\delta_i P_{ij}}



%% begin run
tic
% X0=X0; %迭代初值
% A=A;%投影系数矩阵
% Y=Y;%Tec 测量值

%
% X=X0;

%
% [Nrays,N_net]=size(A);
%
% for j=1:iterations  %没误差10次迭代，有误差3次迭代
%     for i=1:Nrays
%         if dot(X,A(i,:))>0
%             base=Y(i)/dot(X,A(i,:));
%             mo=norm(A(i,:));
%             for l=1:N_net
%                 X(l)=X(l)*(base^(relax*A(i,l)'/mo));
%             end
%         end
%     end
%     fprintf('Iterations(迭代次数) %d\n', j);
% end


iterations=3;  %迭代次数
%
relax=1;  %松弛因子
X= mart(Y,A,X0,relax,iterations);
save('X.mat','X')

toc

%% plot
Nx=52;
Ny=10;
Nz=40;

Ne_true_net=reshape(Ne_true,40,52,10);
Ne_true_2d=reshape(Ne_true_net(:,:,1),40,52);


v=0:.12e11:1.56e11;
Azi=0:0.0288*82.4891:1.5*82.4891-0.0288*82.4891;
figure
contourf(Azi(9:44)-Azi(9),205:5:400,Ne_true_2d(:,9:44),v);
h=colorbar;
set(get(h,'title'),'string');
axis([0 80 205 400])
xlabel('Azimuth Direcction/km');
ylabel('Altitude/km');


X_net=reshape(X,40,52,10);   %反演得到的电子密度分布
X_2d=reshape(X_net(:,:,1),40,52);
figure
contourf(Azi(9:44)-Azi(9),205:5:400,X_2d(:,9:44),v);
h=colorbar;
set(get(h,'title'),'string');
axis([0 80 205 400])
xlabel('Azimuth Direcction/km');
ylabel('Altitude/km');



% [x,y,z] = meshgrid(5:2:109-2,5:20:205-20,100:10:500-10);
[x,y,z] = meshgrid(5:2:109-2,100:10:500-10,5:20:205-20);


h=slice(x,y,z,X_net,[],[],[5:40:205-20]);

for i=1:length(h)
    h(i).FaceColor = 'interp';
    h(i).EdgeColor = 'none';
    h(i).DiffuseStrength = 0.8;
end
colorbar
caxis([0 1.56e11]);


% save('Ne_true.mat','Ne_true');
% save('X.mat','X')
