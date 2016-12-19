function [xk] = RANDOM_KACZMARZ(A, b, x0, it)
% RANDOM_KACMARZ finding a approximate solution to a linear system using the
% radonmized kazcmarz method
% 
% Input
%   A:      System matrix   (Matrix, dim mxn)
%   b:      Right hand side (Vector, dim m)
%   x0:     Initial guess   (Vector, dim n)
%   it:      Number of iterations (Scalar)
%
% Output
%   xk:     The k'th iteration
%
% Authors: (Date 26.04.2016)
%   Sebastian Balle,         s144243
%   Tor Schneider Dengsøe,   s144236
%   Mikkel Paltorp Schmitt,  s144251

xk = x0;
n = length(x0);

probs = sum(abs(A).^2,2)/norm(A, 'fro')^2;
cdf = cumsum(probs);

for i = 1:it
    j = mod(sum(cdf <= rand), n)+1;
    Ai = A(j,:);
    xk1 = xk + ((b(j)-Ai*xk)/norm(Ai,2)^2*Ai)';
    xk = xk1;
end
end