function [A, B, C, D]=SIM_PCA(y,u,i,n)
%利用子空间方法（SIMPCA）对输入输出数据进行辨识
%参考论文 A new subspace identification approach based on principal component analysis
%   y是输出数据
%   u是输入数据
%   i是历史数据的数量，也就是分块Hankel矩阵行数的1/2
%   n是系统阶数
%   A,B,C,D是辨识得到的状态空间矩阵

%   TODO：增加阶次辨识方法

addpath(genpath('../'));

%   判断输入参数
if (nargin < 4);
    error('输入参数不足');
end

[m, num_u] = size(u);
if(m > num_u)
    u = u';
    [m, num_u] = size(u);
end

[l, num_y] = size(y);
if(l > num_y)
    y = y';
    [l, num_y] = size(y);
end

if(num_u ~= num_y)
    error('输入输出向量长度必须相等');
end

num = num_u;
j = num - 2 * i + 1;

%   计算Hankel矩阵
Y = blkhank(y, 2 * i, j);
U = blkhank(u, 2 * i, j);
Yp = Y(1:i * l, :);
Yf = Y(i * l + 1:2 * i * l, :);
Up = U(1:i * m, :);
Uf = U(i * m + 1:2 * i * m, :);
Zp = [Yp; Up];  %注意与基本的SIM不同
Zf = [Yf; Uf];

%利用SVD进行PCA分解（参考论文An Extented Closed-loop Subspace Identification Method for Error-in-variables System 式（27））
[U, S, V] = svd((Zf * (Zp)') / j);


P_resid = U(:, (m * i) + n + 1:end); %提取残余项的载荷矩阵
P_resid_y = P_resid(1:l * i, 1:l*i-n);
P_resid_u = P_resid(l * i + 1:(l + m) * i, 1:l*i-n);
gammaf_per = P_resid_y'; %式（28）中，M取单位阵
gammaf = orthcomp(gammaf_per'); %计算正交补
%gammaf = null(gammaf_per); %计算正交补的另一种方法，数值结果不同，最终A,B,C,D阵也不同，但系统特性基本一样

%   计算A和C
C = gammaf(1:l, 1:n);
A = pinv(gammaf(1:l * (i - 1), 1:n)) * gammaf(l + 1:l * i, 1:n);

%   准备计算B和D，参考式（32）-（36）
phi = -P_resid_y';
Phi = P_resid_u';
Lhs = zeros((l * i - n) * i, l * i);
Rhs = zeros((l * i - n) * i, m);
for k = 1 : i
    Rhs((k - 1) * (l * i - n) + 1:k * (l * i - n), :) = Phi(:, (k - 1) * m + 1:k * m);
    Lhs((k - 1) * (l * i - n) + 1:k * (l * i - n), 1:(i - k + 1) * l) = phi(:,(k - 1) * l + 1:l * i);
end

%   求解最小二乘
Hf = Lhs \ Rhs;
sol = pinv([eye(l) zeros(l, n); zeros(l * (i - 1), l) gammaf(1:l * (i - 1), 1:n)]) * Hf; %式（40）

%   得到B和D矩阵
D = sol(1:l, :);
B = sol(l + 1:l + n, :);
