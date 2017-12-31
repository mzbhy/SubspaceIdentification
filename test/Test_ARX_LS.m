%测试基于ARX模型的最小二乘辨识算法
%   VarName1是从SD卡中导入的MTi原始数据（辨识用），为float型单列向量
%   VarName2是从SD卡中导入的AD7606原始数据（辨识用），为单列向量
%   VarName3是从SD卡中导入的MTi原始数据（验证用），为float型单列向量
%   VarName4是从SD卡中导入的AD7606原始数据（验证用），为单列向量
%   M是最小二乘法得到的系统模型

addpath(genpath('../'));
clearvars -EXCEPT VarName1 VarName2 VarName3 VarName4

%%%%采样信息设置%%%%%
T = 0.02;
Fs = 1 / T;              % 采样频率

%%%%%删除NAN数据%%%%%
NaNlocation = find(isnan(VarName1));
VarName1(NaNlocation)=[];
VarName2(NaNlocation)=[];

%%%%%获取输入输出数据%%%%%
[beta, diff_MTi, beta_filter] = Get_JY01_AD7606_MTi(VarName1, VarName2, T);

%%%%%最小二乘辨识%%%%%
data_length = min(length(diff_MTi), length(beta_filter));
M = ARX_LS_SYSID(diff_MTi(1:data_length), beta_filter(1:data_length), Fs, round(0.8*length(diff_MTi)));

%%%%%如果需要其他数据集进行验证，取消下面的注释%%%%%
%%%%%删除NAN数据%%%%%
NaNlocation = find(isnan(VarName3));
VarName3(NaNlocation)=[];
VarName4(NaNlocation)=[];

%%%%%获取输入输出数据%%%%%
[beta_verify, diff_MTi_verify, beta_filter_verify] = Get_JY01_AD7606_MTi(VarName3, VarName4, T);

%%%%%测试输出%%%%%
YP = ARX_LS_SYSID_Verify(diff_MTi_verify, beta_filter_verify(1:length(diff_MTi_verify)), Fs, M);
plot(YP,'r');
hold on
plot(beta_filter_verify,'b');
legend('辨识输出','真实输出');
title('辨识模型的输出');

