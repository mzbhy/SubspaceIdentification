function [] = Id_PZ_Plot(sys)
%绘制零极点图示例
%   sys是待绘制系统

figure(1)
pzmap(sys,'b');
grid on