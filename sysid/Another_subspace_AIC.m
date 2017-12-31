function [AIC] = Another_subspace_AIC(u, y, ss, n)
%另一种AIC定义，参考《Subspace identification using the integration of MOSEP and N4SID methods applied to the Shell benchmark of a distillation column》利用最小二乘法，采用ARX模型，对系统进行辨识。
%   u是输入向量
%   y是实际输出向量
%   ss是阶次n下得到的辨识系统
%   n是系统阶次

nu = length(u);
y_hat = dlsim(ss.A, ss.B, ss.C, ss.D, u);
e = y - y_hat;

AIC = nu * log(var(e)) + 4 * nu * (n^2 + 2 * n) / (nu - (n^2 + 2 * n + 1));

end