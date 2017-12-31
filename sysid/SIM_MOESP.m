function [A, B, C, D] = SIM_MOESP(y, u, i, n)
%利用子空间方法(Multivariable Output Error State Space)对输入输出数据进行辨识
%参考Subspace identification for System Identification 159页。
%   y是输出数据
%   u是输入数据
%   i是分块Hankel矩阵行数
%   n是系统阶数，可以不指定
%   A,B,C,D是辨识得到的状态空间矩阵

addpath(genpath('../'));

%   判断输入参数，如果n未制定，准备进行阶次辨识
if (nargin < 4);
    n = [];
end

[num_u, m]=size(u);
if(m > num_u)
    u = u';
    [num_u, m]=size(u);
end

[num_y, l]=size(y);
if(l > num_y)
    y = y';
    [num_y, l]=size(y);
end

if num_u~=num_y
    error('输入输出向量长度必须相等')
end

%   计算Hankel矩阵
Y = blkhank(y, i, num_y - i + 1);
U = blkhank(u, i, num_y - i + 1);

%   LQ分解
L = triu(qr([U;Y]'))';
L = L(1:i * (m + l), 1:i * (m + l));

%   SVD
L22 = L(i*l+1:end,i*l+1:end);
[U1,S1]=svd(L22);
ss = diag(S1);

%   阶次辨识
if(isempty(n))
    x = 1:1:l * i;
    figure();
    subplot;
    bar(x, ss);
    axis([0, length(ss) + 1, 10^(floor(log10(min(ss)))), 10^(ceil(log10(max(ss))))]);
    title('奇异值');
    xlabel('阶次');
    n = input('系统阶次应为？');
end


% C and A
Ok = U1(:, 1:n) * diag(sqrt(ss(1:n)));
C = Ok(1:l, :);
A = Ok(1:l * (i - 1), :) \ Ok(l + 1:i * l, :);

% B and D 式(6.43)、(6.44)
U2 = U1(:, n + 1:end)';
L11 = L(1:i * m, 1:i * m);
L21 = L(i * m + 1:end,1:i * m);
M1 = U2 * L21 / L11;
M = zeros((l * i - n) * i, m);
L = zeros((l * i - n) * i, l + n);
for k = 1:i
    M((k - 1) * (l * i - n) + 1:k * (l * i - n), :) = M1(:, (k - 1) * m + 1:k * m);
    L((k - 1) * (l * i - n) + 1:k * (l * i - n), :) = [U2(:, (k - 1) * l + 1:k * l) U2(:, k * l + 1:end) * Ok(1:end - k * l, :)];
end

DB = L \ M;
D = DB(1 : l, :);
B = DB(l + 1:end, :);
end
