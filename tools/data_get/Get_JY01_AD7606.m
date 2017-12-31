function [beta, beta_filter] = Get_JY01_AD7606(VarName1, T)
%用于得到利用AD7606采集的JY01数据，并进行初步处理。
%   VarName1是从SD卡中导入的AD7606原始数据，为单列向量
%   T是采样时间
%   beta是AD数据直接转换得到的角加速度数据
%   beta_filter是对beta进行滤波后得到的数据

addpath(genpath('../'));

%%%%采样信息设置%%%%%
Fs = 1 / T;              % 采样频率
L = length(VarName1);     % 信号长度
time = T * (0:1:L-1);


%%%%%AD原始数据（JY01）处理%%%%%
beta = JY01_ADValue_Process(VarName1, 5);
figure;
plot(time, beta,'r');
title('角加速度信号');

%%%%%AD采用巴特沃斯滤波器%%%%%
Rp = 1;As = 30;
fp = 10;fs = 20;
[beta_filter, H] = Butterworth_Filter(beta, Fs, fp, fs, Rp, As);

% figure
% plot(W, abs(H))
% title('滤波器频率特性');

figure
plot(time, beta_filter);
title('JY01滤波结果');

%%%%%计算JY01的fft%%%%%
Cal_FFT(beta, Fs);
