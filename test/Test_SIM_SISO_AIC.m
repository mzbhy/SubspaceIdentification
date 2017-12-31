%SISO仿真数据，测试各类型的子空间算法，并绘制伯德图
%   SIM_Deterministic是Subspace identification for linear system中的基本子空间算法
%   SIM_MOESP是Subspace identification for System Identification中的Multivariable Output Error State Space算法
%   SIM_PCA是A new subspace identification approach based on principal component analysis中的基于PCA的子空间算法

clear;
addpath(genpath('../'));

N = 1000;
u = randn(N, 1);
u_test =u + 0 * randn(N, 1);
a = [0.6 0.4;-0.4 0.6];
b = [0;1];
c = [1 0.5];
d = [0];
y = dlsim(a, b, c, d, u_test);
y_test = y + 0 * randn(N, 1);
k = 10;

% [A_SIM_D, B_SIM_D, C_SIM_D, D_SIM_D] = SIM_Deterministic(y_test, u_test, 100, 2);
% [A_SIM_MOESP, B_SIM_MOESP, C_SIM_MOESP, D_SIM_MOESP] = SIM_MOESP(y_test, u_test, 200, 2);
% [A_SIM_PCA, B_SIM_PCA, C_SIM_PCA, D_SIM_PCA] = SIM_PCA(y_test, u_test, 100, 2);
% sys_real = ss(a,b,c,d);
% sys_sim_d = ss(A_SIM_D, B_SIM_D, C_SIM_D, D_SIM_D);
% sys_sim_mosep = ss(A_SIM_MOESP, B_SIM_MOESP, C_SIM_MOESP, D_SIM_MOESP);
% sys_sim_pca = ss(A_SIM_PCA, B_SIM_PCA, C_SIM_PCA, D_SIM_PCA);



[A_SIM_MOESP, B_SIM_MOESP, C_SIM_MOESP, D_SIM_MOESP] = SIM_MOESP(y_test, u_test, 20, 1);
[A_SIM_MOESP2, B_SIM_MOESP2, C_SIM_MOESP2, D_SIM_MOESP2] = SIM_MOESP(y_test, u_test, 20, 2);
[A_SIM_MOESP3, B_SIM_MOESP3, C_SIM_MOESP3, D_SIM_MOESP3] = SIM_MOESP(y_test, u_test, 20, 3);
[A_SIM_MOESP4, B_SIM_MOESP4, C_SIM_MOESP4, D_SIM_MOESP4] = SIM_MOESP(y_test, u_test, 20, 4);
[A_SIM_MOESP5, B_SIM_MOESP5, C_SIM_MOESP5, D_SIM_MOESP5] = SIM_MOESP(y_test, u_test, 20, 5);
[A_SIM_MOESP6, B_SIM_MOESP6, C_SIM_MOESP6, D_SIM_MOESP6] = SIM_MOESP(y_test, u_test, 20, 6);
[A_SIM_MOESP7, B_SIM_MOESP7, C_SIM_MOESP7, D_SIM_MOESP7] = SIM_MOESP(y_test, u_test, 20, 7);
[A_SIM_MOESP8, B_SIM_MOESP8, C_SIM_MOESP8, D_SIM_MOESP8] = SIM_MOESP(y_test, u_test, 20, 8);
[A_SIM_MOESP9, B_SIM_MOESP9, C_SIM_MOESP9, D_SIM_MOESP9] = SIM_MOESP(y_test, u_test, 20, 9);
[A_SIM_MOESP10, B_SIM_MOESP10, C_SIM_MOESP10, D_SIM_MOESP10] = SIM_MOESP(y_test, u_test, 20, 10);

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
AIC(1) = subspace_AIC(A_SIM_MOESP, B_SIM_MOESP, C_SIM_MOESP, D_SIM_MOESP, y_test, u_test, 20, 1);
AIC(2) = subspace_AIC(A_SIM_MOESP2, B_SIM_MOESP2, C_SIM_MOESP2, D_SIM_MOESP2, y_test, u_test, 20, 2);
AIC(3) = subspace_AIC(A_SIM_MOESP3, B_SIM_MOESP3, C_SIM_MOESP3, D_SIM_MOESP3, y_test, u_test, 20, 3);
AIC(4) = subspace_AIC(A_SIM_MOESP4, B_SIM_MOESP4, C_SIM_MOESP4, D_SIM_MOESP4, y_test, u_test, 20, 4);
AIC(5) = subspace_AIC(A_SIM_MOESP5, B_SIM_MOESP5, C_SIM_MOESP5, D_SIM_MOESP5, y_test, u_test, 20, 5);
AIC(6) = subspace_AIC(A_SIM_MOESP6, B_SIM_MOESP6, C_SIM_MOESP6, D_SIM_MOESP6, y_test, u_test, 20, 6);
AIC(7) = subspace_AIC(A_SIM_MOESP7, B_SIM_MOESP7, C_SIM_MOESP7, D_SIM_MOESP7, y_test, u_test, 20, 7);
AIC(8) = subspace_AIC(A_SIM_MOESP8, B_SIM_MOESP8, C_SIM_MOESP8, D_SIM_MOESP8, y_test, u_test, 20, 8);
AIC(9) = subspace_AIC(A_SIM_MOESP9, B_SIM_MOESP9, C_SIM_MOESP9, D_SIM_MOESP9, y_test, u_test, 20, 9);
AIC(10) = subspace_AIC(A_SIM_MOESP10, B_SIM_MOESP10, C_SIM_MOESP10, D_SIM_MOESP10, y_test, u_test, 20, 10);

figure()
plot(AIC);
% 
% anotherAIC = zeros(1,10);
% anotherAIC(1) = Another_subspace_AIC(y_test, u_test, M_SIM_MOESP, 1);
% anotherAIC(2) = Another_subspace_AIC(y_test, u_test, M_SIM_MOESP2, 2);
% anotherAIC(3) = Another_subspace_AIC(y_test, u_test, M_SIM_MOESP3, 3);
% anotherAIC(4) = Another_subspace_AIC(y_test, u_test, M_SIM_MOESP4, 4);
% anotherAIC(5) = Another_subspace_AIC(y_test, u_test, M_SIM_MOESP5, 5);
% anotherAIC(6) = Another_subspace_AIC(y_test, u_test, M_SIM_MOESP6, 6);
% anotherAIC(7) = Another_subspace_AIC(y_test, u_test, M_SIM_MOESP7, 7);
% anotherAIC(8) = Another_subspace_AIC(y_test, u_test, M_SIM_MOESP8, 8);
% anotherAIC(9) = Another_subspace_AIC(y_test, u_test, M_SIM_MOESP9, 9);
% anotherAIC(10) = Another_subspace_AIC(y_test, u_test, M_SIM_MOESP10, 10);
% 
% figure()
% plot(anotherAIC);

