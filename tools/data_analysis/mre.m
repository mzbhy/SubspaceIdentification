function [mre_value] = mre(y_test, y)
%计算平均相对误差
%   y_test是测量值向量
%   y是真值向量
%   mre_value是测量值与真值之间的平均相对误差

e = y_test - y;
mre_value = sqrt(sum(e.^2) / sum(y.^2));

end

