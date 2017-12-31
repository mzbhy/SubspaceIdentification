%MIMO仿真数据，测试各类型的子空间算法，并绘制伯德图
%   SIM_Deterministic是Subspace identification for linear system中的基本子空间算法
%   SIM_MOESP是Subspace identification for System Identification中的Multivariable Output Error State Space算法
%   SIM_PCA是A new subspace identification approach based on principal component analysis中的基于PCA的子空间算法

clear;
addpath(genpath('../'));

N = 10000;
u = randn(N, 2);
u_test =u + 0.1 * randn(N, 2);
% a = [0.603 0.603 0 0;-0.603 0.603 0 0;0 0 -0.603 -0.603;0 0 0.603 -0.603];
% b = [1.1650,-0.6965;0.6268 1.6961;0.0751,0.0591;0.3516 1.7971];
% c = [0.2641,-1.4462,1.2460,0.5774];
% d = [-0.1356,-1.2704];
a = [0.603 0.603 0 0;-0.603 0.603 0 0;0 0 -0.603 -0.603;0 0 0.603 -0.603];
b = [1.1650,-0.6965;0.6268 1.6961;0.0751,0.0591;0.3516 1.7971];
c = [0.2641,-1.4462,1.2460,0.5774;0.8717,-0.7012,-0.6390,-0.3600];
d = [-0.1356,-1.2704;-1.3493,0.9846];
y = dlsim(a, b, c, d, u_test);
y_test = y + 0.1 * randn(N, 2);
k = 10;

[A_SIM_D, B_SIM_D, C_SIM_D, D_SIM_D] = SIM_Deterministic(y_test, u_test, 100);
% [A_SIM_MOESP, B_SIM_MOESP, C_SIM_MOESP, D_SIM_MOESP] = SIM_MOESP(y_test, u_test, 200, 4);
[A_SIM_PCA, B_SIM_PCA, C_SIM_PCA, D_SIM_PCA] = SIM_PCA(y_test, u_test, 100, 4);
sys_real = ss(a,b,c,d);
sys_sim_d = ss(A_SIM_D, B_SIM_D, C_SIM_D, D_SIM_D);
% sys_sim_mosep = ss(A_SIM_MOESP, B_SIM_MOESP, C_SIM_MOESP, D_SIM_MOESP);
sys_sim_pca = ss(A_SIM_PCA, B_SIM_PCA, C_SIM_PCA, D_SIM_PCA);



figure(1)
h = bodeplot(sys_real,'r',sys_sim_d,'*g',sys_sim_pca,'y');%,sys_sim_mosep,'b',sys_sim_pca,'y');
p = getoptions(h); 
p.PhaseMatching = 'on'; 
setoptions(h,p); % Update the Bode plot
legend('真值', '子空间辨识模型', 'SIM-PCA');%, 'MOESP', 'SIM-PCA');

% figure(2)
% pzmap(sys_real,'r');
% hold on
% pzmap(sys_sim_d,'b');
% legend('真值', '子空间辨识模型');

y_SIM_D = dlsim(A_SIM_D, B_SIM_D, C_SIM_D, D_SIM_D, u_test);
y_SIM_PCA = dlsim(A_SIM_PCA, B_SIM_PCA, C_SIM_PCA, D_SIM_PCA, u_test);
figure(3)
subplot(2,1,1); 
plot(y_test(:,1), 'r','linewidth',3);
hold on
plot(y_SIM_D(:,1), 'b','linewidth',2);
subplot(2,1,2); 
plot(y_test(:,2), 'r','linewidth',3);
hold on
plot(y_SIM_D(:,2), 'b','linewidth',2);
legend('真值', '子空间辨识模型');
figure(4)
subplot(2,1,1); 
plot(y_test(:,1), 'r','linewidth',3);
hold on
plot(y_SIM_PCA(:,1), 'b','linewidth',2);
subplot(2,1,2); 
plot(y_test(:,2), 'r','linewidth',3);
hold on
plot(y_SIM_PCA(:,2), 'b','linewidth',2);
legend('真值', '子空间辨识模型');
