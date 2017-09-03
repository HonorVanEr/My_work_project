function [ X ] = mart(Y,A,X0, relax, iterations) %#ok<INUSD>

%mart Performs multiplicitive algebraic reconstruction (MART)
%
% Usage: [X] = mart(Y,A,X0,relax,iterations)
%

% Inputs:
%   -Y: Array containing set of projection data. TEC data Nrays x 1
%   -A: Projection matrix for Y. from matrix_3D Nrays x N_net
%   -X0: Initial guess for Ne.������ֵ IRI model N_net x 1
%   -relax:  relaxation parameter default 1.0
%   -iterations:   The number of iterations to perform during
%   reconstruction.��������
%
% Ouputs:
%   -X: Ne N_net x 1
%
% Notes: must assume positivity of image--which makes sense for ionosphere tomography, thus x0 must be > 0
%     x_j^{k+1} = x_j^k * (y_i / (\sum_j a_{ij}x_j^k))^{\gamma\delta_i P_{ij}}
%
%     Implementation: see "The Multiplicative Algebraic Reconstruction Technique Solves the Geometric Programming Problem"
%         by Charles Byrne, October 23, 2007

X=X0;
% iterations=3;  %��������
% relax=1;  %�ɳ�����

[Nrays,N_net]=size(A);

for j=1:iterations  %û���10�ε����������3�ε���
    for i=1:length(A(:,1))
        if dot(X,A(i,:))>0
            base=Y(i)/dot(X,A(i,:));
            mo=norm(A(i,:));
            for l=1:N_net
                X(l)=X(l)*(base^(relax*A(i,l)'/mo));
            end
        end
    end
    fprintf('Iterations(��������) %d\n', j);
end

end

