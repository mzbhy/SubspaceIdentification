function [u_test, y_test] = Monte_Carlo_Simulation(a, b, c, d, N, err_u, err_y)
%产生测试数据，对辨识算法进行测试
%   a,b,c,d是测试系统的系统矩阵
%   N是产生数据长度
%   err_u是输入信号的噪声标准差
%   err_y是输出信号的噪声标准差
%   u_test是产生的输入信号
%   y_test是产生的输出信号

%   判断输入参数
if (nargin < 5);
    error('输入参数不足');
end
if (isequal(nargin, 5))
    err_u = 0.1;
    err_y = 0.1;
elseif (isequal(nargin, 6))
    err_y = 0.1;
end

[n_u, size_u] = size(b);
[size_y, n_y] = size(c);
if(n_u ~= n_y)
    error('系统矩阵维数不匹配');
end


u = randn(N, size_u);
u_test =u + err_u * randn(N, size_u);
y = dlsim(a, b, c, d, u_test);
y_test = y + err_y * randn(N, size_y);

