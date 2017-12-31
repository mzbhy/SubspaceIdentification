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

[A_SIM_MOESP, B_SIM_MOESP, C_SIM_MOESP, D_SIM_MOESP] = SIM_MOESP(beta_filter(1:data_length), diff_MTi(1:data_length), 40, 1);
[A_SIM_MOESP2, B_SIM_MOESP2, C_SIM_MOESP2, D_SIM_MOESP2] = SIM_MOESP(beta_filter(1:data_length), diff_MTi(1:data_length), 40, 2);
[A_SIM_MOESP3, B_SIM_MOESP3, C_SIM_MOESP3, D_SIM_MOESP3] = SIM_MOESP(beta_filter(1:data_length), diff_MTi(1:data_length), 40, 3);
[A_SIM_MOESP4, B_SIM_MOESP4, C_SIM_MOESP4, D_SIM_MOESP4] = SIM_MOESP(beta_filter(1:data_length), diff_MTi(1:data_length), 40, 4);
[A_SIM_MOESP5, B_SIM_MOESP5, C_SIM_MOESP5, D_SIM_MOESP5] = SIM_MOESP(beta_filter(1:data_length), diff_MTi(1:data_length), 40, 5);
[A_SIM_MOESP6, B_SIM_MOESP6, C_SIM_MOESP6, D_SIM_MOESP6] = SIM_MOESP(beta_filter(1:data_length), diff_MTi(1:data_length), 40, 6);
[A_SIM_MOESP7, B_SIM_MOESP7, C_SIM_MOESP7, D_SIM_MOESP7] = SIM_MOESP(beta_filter(1:data_length), diff_MTi(1:data_length), 40, 7);
[A_SIM_MOESP8, B_SIM_MOESP8, C_SIM_MOESP8, D_SIM_MOESP8] = SIM_MOESP(beta_filter(1:data_length), diff_MTi(1:data_length), 40, 8);
[A_SIM_MOESP9, B_SIM_MOESP9, C_SIM_MOESP9, D_SIM_MOESP9] = SIM_MOESP(beta_filter(1:data_length), diff_MTi(1:data_length), 40, 9);
[A_SIM_MOESP10, B_SIM_MOESP10, C_SIM_MOESP10, D_SIM_MOESP10] = SIM_MOESP(beta_filter(1:data_length), diff_MTi(1:data_length), 40, 10);

M_SIM_MOESP = ss(A_SIM_MOESP, B_SIM_MOESP, C_SIM_MOESP, D_SIM_MOESP);
M_SIM_MOESP2 = ss(A_SIM_MOESP2, B_SIM_MOESP2, C_SIM_MOESP2, D_SIM_MOESP2);
M_SIM_MOESP3 = ss(A_SIM_MOESP3, B_SIM_MOESP3, C_SIM_MOESP3, D_SIM_MOESP3);
M_SIM_MOESP4 = ss(A_SIM_MOESP4, B_SIM_MOESP4, C_SIM_MOESP4, D_SIM_MOESP4);
M_SIM_MOESP5 = ss(A_SIM_MOESP5, B_SIM_MOESP5, C_SIM_MOESP5, D_SIM_MOESP5);
M_SIM_MOESP6 = ss(A_SIM_MOESP6, B_SIM_MOESP6, C_SIM_MOESP6, D_SIM_MOESP6);
M_SIM_MOESP7 = ss(A_SIM_MOESP7, B_SIM_MOESP7, C_SIM_MOESP7, D_SIM_MOESP7);
M_SIM_MOESP8 = ss(A_SIM_MOESP8, B_SIM_MOESP8, C_SIM_MOESP8, D_SIM_MOESP8);
M_SIM_MOESP9 = ss(A_SIM_MOESP9, B_SIM_MOESP9, C_SIM_MOESP9, D_SIM_MOESP9);
M_SIM_MOESP10 = ss(A_SIM_MOESP10, B_SIM_MOESP10, C_SIM_MOESP10, D_SIM_MOESP10);

AIC = zeros(1,10);
AIC(1) = subspace_AIC(A_SIM_MOESP, B_SIM_MOESP, C_SIM_MOESP, D_SIM_MOESP, beta_filter(1:data_length), diff_MTi(1:data_length), 20, 1);
AIC(2) = subspace_AIC(A_SIM_MOESP2, B_SIM_MOESP2, C_SIM_MOESP2, D_SIM_MOESP2, beta_filter(1:data_length), diff_MTi(1:data_length), 20, 2);
AIC(3) = subspace_AIC(A_SIM_MOESP3, B_SIM_MOESP3, C_SIM_MOESP3, D_SIM_MOESP3, beta_filter(1:data_length), diff_MTi(1:data_length), 20, 3);
AIC(4) = subspace_AIC(A_SIM_MOESP4, B_SIM_MOESP4, C_SIM_MOESP4, D_SIM_MOESP4, beta_filter(1:data_length), diff_MTi(1:data_length), 20, 4);
AIC(5) = subspace_AIC(A_SIM_MOESP5, B_SIM_MOESP5, C_SIM_MOESP5, D_SIM_MOESP5, beta_filter(1:data_length), diff_MTi(1:data_length), 20, 5);
AIC(6) = subspace_AIC(A_SIM_MOESP6, B_SIM_MOESP6, C_SIM_MOESP6, D_SIM_MOESP6, beta_filter(1:data_length), diff_MTi(1:data_length), 20, 6);
AIC(7) = subspace_AIC(A_SIM_MOESP7, B_SIM_MOESP7, C_SIM_MOESP7, D_SIM_MOESP7, beta_filter(1:data_length), diff_MTi(1:data_length), 20, 7);
AIC(8) = subspace_AIC(A_SIM_MOESP8, B_SIM_MOESP8, C_SIM_MOESP8, D_SIM_MOESP8, beta_filter(1:data_length), diff_MTi(1:data_length), 20, 8);
AIC(9) = subspace_AIC(A_SIM_MOESP9, B_SIM_MOESP9, C_SIM_MOESP9, D_SIM_MOESP9, beta_filter(1:data_length), diff_MTi(1:data_length), 20, 9);
AIC(10) = subspace_AIC(A_SIM_MOESP10, B_SIM_MOESP10, C_SIM_MOESP10, D_SIM_MOESP10, beta_filter(1:data_length), diff_MTi(1:data_length), 20, 10);

figure()
plot(AIC);

anotherAIC = zeros(1,10);
anotherAIC(1) = Another_subspace_AIC(diff_MTi(1:data_length), beta_filter(1:data_length), M_SIM_MOESP, 1);
anotherAIC(2) = Another_subspace_AIC(diff_MTi(1:data_length), beta_filter(1:data_length), M_SIM_MOESP2, 2);
anotherAIC(3) = Another_subspace_AIC(diff_MTi(1:data_length), beta_filter(1:data_length), M_SIM_MOESP3, 3);
anotherAIC(4) = Another_subspace_AIC(diff_MTi(1:data_length), beta_filter(1:data_length), M_SIM_MOESP4, 4);
anotherAIC(5) = Another_subspace_AIC(diff_MTi(1:data_length), beta_filter(1:data_length), M_SIM_MOESP5, 5);
anotherAIC(6) = Another_subspace_AIC(diff_MTi(1:data_length), beta_filter(1:data_length), M_SIM_MOESP6, 6);
anotherAIC(7) = Another_subspace_AIC(diff_MTi(1:data_length), beta_filter(1:data_length), M_SIM_MOESP7, 7);
anotherAIC(8) = Another_subspace_AIC(diff_MTi(1:data_length), beta_filter(1:data_length), M_SIM_MOESP8, 8);
anotherAIC(9) = Another_subspace_AIC(diff_MTi(1:data_length), beta_filter(1:data_length), M_SIM_MOESP9, 9);
anotherAIC(10) = Another_subspace_AIC(diff_MTi(1:data_length), beta_filter(1:data_length), M_SIM_MOESP10, 10);

figure()
plot(anotherAIC);
