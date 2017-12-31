function [rmse_value] = rmse(y_test, y)
%计算均方根误差
%   y_test是测量值向量
%   y是真值向量
%   rmse_value是测量值与真值之间的均方值误差

rmse_value = sqrt(mean((y_test - y).^2));

end

