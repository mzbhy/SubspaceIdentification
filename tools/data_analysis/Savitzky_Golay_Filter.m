function y = Savitzky_Golay_Filter(u, N, F)
%Savitzky-Golay Filter.
%   u是原始数据，为单行（或单列）向量
%   N是滤波器阶数
%   F是滤波器窗口宽度，必须为奇数
%   y是滤波后数据，为单行（或单列）向量，长度比u少(F+1)/2

[~, g] = sgolay(N, F);   % 计算S-G滤波器参数

HalfWin  = ((F + 1) / 2) - 1; %半窗口宽度
y = zeros(1, length(u) - (F+1)/2);
for n = (F + 1) / 2 : length(u) - (F + 1) / 2
  y(n) = dot(g(:, 1), u(n - HalfWin : n + HalfWin));
end
%plot([u(1 : length(y)), y']);