function [beta, diff_MTi, beta_filter, VarName1] = Get_JY01_AD7606_MTi(VarName1, VarName2, T, filter_method)
%用于得到利用AD7606采集的JY01数据，以及利用MTi惯性导航模块得到的陀螺仪数据，进行初步处理。
%   VarName1是从SD卡中导入的MTi原始数据，为float型单列向量
%   VarName2是从SD卡中导入的AD7606原始数据，为单列向量
%   T是采样时间
%   如果filter_method是对JY01信号进行滤波的方法，可选项有：'butterworth','S-G','kalman','FIR'，默认值为'butterworth'
%   beta是AD数据直接转换得到的角加速度数据
%   beta_filter是对beta进行滤波后得到的数据
%   diff_MTi是MTi数据微分得到的角加速度数据
%   VarName1是MTi原始数据

addpath(genpath('../'));

%   判断输入参数，如果filter_method未制定，默认为巴特沃斯低通滤波器
if (nargin < 4);
    filter_method = 'butterworth';
end

%%%%采样信息设置%%%%%
Fs = 1 / T;              % 采样频率
L = length(VarName1);     % 信号长度
time = T * (0:1:L-1);

%%%%%Mti原始数据%%%%%
plot(time, VarName1)
title('Mti原始数据');

% %%%%%MTi采用巴特沃斯滤波器%%%%%
% Rp = 1;As = 30;
% fp = 10;fs = 20;
% [MTi_filter, H_MTi] = Butterworth_Filter(VarName1, Fs, fp, fs, Rp, As);

%%%%%Mti数据处理%%%%%
diff_MTi = MTi_Process(VarName1, 'center')  * 180 / pi / T;
figure
plot(time(1:length(diff_MTi)), diff_MTi, 'b');


%%%%%AD原始数据（JY01）处理%%%%%
beta = JY01_ADValue_Process(VarName2, 5);
hold on;
plot(time, beta,'r');
legend('MTi角加速度信号','JY01角加速度信号');
title('两种角加速度信号');

%%%%%滤波器设置%%%%%
switch(filter_method)
    case 'butterworth'
        Rp = 1;As = 30;
        fp = 5;fs = 20;
        [beta_filter, H] = Butterworth_Filter(beta, Fs, fp, fs, Rp, As);
        % figure
        % plot(W, abs(H))
        % title('滤波器频率特性');
    case 'S-G'
        beta_filter = Savitzky_Golay_Filter(beta, 3, 15);
    case 'kalman'
    case 'FIR'
        Rp = 1;As = 30;
        fp = 5;fs = 20;
        beta_filter = FIR_Filter(beta, Fs, fp, fs, As);
end

figure
plot(time(1:length(beta_filter)), beta_filter);
title('JY01滤波结果');

%%%%%计算JY01的fft%%%%%
Cal_FFT(beta, Fs);

%%%%%计算误差%%%%%
err = diff_MTi - beta(1:length(diff_MTi));
figure
plot(time(1:length(diff_MTi)), err);
title('两种角加速度信号的误差');
