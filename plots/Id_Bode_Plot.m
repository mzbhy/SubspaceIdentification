function [] = Id_Bode_Plot(sys)
%绘制Bode图示例
%   sys是待绘制系统

figure(1)
h = bodeplot(sys, 'b');
p = getoptions(h); 
p.PhaseMatching = 'on';
p.FreqUnits = 'Hz';
setoptions(h,p);