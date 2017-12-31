function N = orthcomp(A)
%求解矩阵的正交补
%   A是待求矩阵
%   N是矩阵A的正交补

A=A';
[R, pivcol] = rref(A, sqrt(eps)); %)将A化成阶梯状
[m, n] = size(A);
r = length(pivcol);
freecol = 1:n;
freecol(pivcol) = [];
N = zeros(n, n-r);
N(freecol, : ) = eye(n-r);
N(pivcol,  : ) = -R(1:r, freecol);