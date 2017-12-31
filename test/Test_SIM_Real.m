 %角加实际数据，测试各类型的子空间算法
%   VarName1是从SD卡中导入的MTi原始数据（辨识用），为float型单列向量
%   VarName2是从SD卡中导入的AD7606原始数据（辨识用），为单列向量
%   VarName3是从SD卡中导入的MTi原始数据（验证用），为float型单列向量
%   VarName4是从SD卡中导入的AD7606原始数据（验证用），为单列向量
%   M是子空间算法得到的系统模型

addpath(genpath('../'));
clearvars -EXCEPT VarName1 VarName2 VarName3 VarName4

%%%%采样信息设置%%%%%
T = 0.02;
Fs = 1 / T;              % 采样频率

%%%%%删除NAN数据%%%%%
NaNlocation = find(isnan(VarName2));
VarName1(NaNlocation)=[];
VarName2(NaNlocation)=[];

%%%%%获取输入输出数据%%%%%
[beta, diff_MTi, beta_filter] = Get_JY01_AD7606_MTi(VarName1, VarName2, T);

%%%%%SIM辨识%%%%%
data_length = min(length(diff_MTi), length(beta_filter));
[A, B, C, D] = SIM_Deterministic(beta_filter(1:data_length), diff_MTi(1:data_length), 20, 4);
%[A_SIM_MOESP, B_SIM_MOESP, C_SIM_MOESP, D_SIM_MOESP] = SIM_MOESP(beta_filter(1:data_length), diff_MTi(1:data_length), 40, 2);
% [A_SIM_PCA, B_SIM_PCA, C_SIM_PCA, D_SIM_PCA] = SIM_PCA(beta_filter(1:data_length), diff_MTi(1:data_length), 60, 2);
M = ss(A, B, C, D);
% M_SIM_MOESP = ss(A_SIM_MOESP, B_SIM_MOESP, C_SIM_MOESP, D_SIM_MOESP);
% M_SIM_PCA = ss(A_SIM_PCA, B_SIM_PCA, C_SIM_PCA, D_SIM_PCA);

% %%%%%如果需要其他数据集进行验证，取消下面的注释%%%%%
% %%%%%删除NAN数据%%%%%
NaNlocation = find(isnan(VarName4));
VarName3(NaNlocation)=[];
VarName4(NaNlocation)=[];

%%%%%获取输入输出数据%%%%%
[beta_verify, diff_MTi_verify, beta_filter_verify] = Get_JY01_AD7606_MTi(VarName3, VarName4, T);

YP = dlsim(A, B, C, D, diff_MTi_verify);
% YP_SIM_MOESP = dlsim(A_SIM_MOESP, B_SIM_MOESP, C_SIM_MOESP, D_SIM_MOESP, diff_MTi_verify);
% YP_SIM_PCA = dlsim(A_SIM_PCA, B_SIM_PCA, C_SIM_PCA, D_SIM_PCA, diff_MTi_verify);
plot(YP,'r');
hold on
% plot(YP_SIM_MOESP,'o-g');
% plot(YP_SIM_PCA,'*-b');
plot(beta_filter_verify,'k');
legend('SIM辨识输出','真实输出');
title('辨识模型的输出');