function X = art_solver(A, b, I)
% ==============================================================
% art_solver Finding a approximate solution to a linear system using the
% kazcmarz method
% 
% Input
%   A:          System matrix        (Matrix, dim mxn)
%   b:          Right hand side      (Vector, dim m)
%   I:          Number of iterations (Scalar)
%
% Output
%   xk:         The k'th iteration
%
% Authors: (Date 26.04.2016)
%   Sebastian Balle,         s144243
%   Tor Schneider Dengsoe,   s144236
%   Mikkel Paltorp Schmitt,  s144251
% ==============================================================

n = size(A);
m = n(1);
n = n(2);
xk = ones(n,1);

for i = 1:I
    j = mod(i, m)+1;
    Ai = A(j,:);
    xk1 = xk + ((b(j)-Ai*xk)/(Ai*Ai')*Ai)';
    xk = xk1;
end
X = xk;
end