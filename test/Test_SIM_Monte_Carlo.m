%SISO仿真数据，测试各类型的子空间算法，并绘制伯德图
%   SIM_Deterministic是Subspace identification for linear system中的基本子空间算法
%   SIM_MOESP是Subspace identification for System Identification中的Multivariable Output Error State Space算法
%   SIM_PCA是A new subspace identification approach based on principal component analysis中的基于PCA的子空间算法

clear;
addpath(genpath('../'));

times = 10;
N = 10000;
err_u = 0.001;
err_y = 0.001;
a = [0.6 0.4;-0.4 0.6];
b = [0;1];
c = [1 0.5];
d = [0];
% a = [0.603 0.603 0 0;-0.603 0.603 0 0;0 0 -0.603 -0.603;0 0 0.603 -0.603];
% b = [1.1650,-0.6965;0.6268 1.6961;0.0751,0.0591;0.3516 1.7971];
% c = [0.2641,-1.4462,1.2460,0.5774;0.8717,-0.7012,-0.6390,-0.3600];
% d = [-0.1356,-1.2704;-1.3493,0.9846];

sys_real = ss(a,b,c,d);

tic
while times > 0
    times = times -1;
    [u_test, y_test] = Monte_Carlo_Simulation(a, b, c, d, N, err_u, err_y);
    [A_SIM_D, B_SIM_D, C_SIM_D, D_SIM_D] = SIM_Deterministic(y_test, u_test, 100, 2);
    
    sys_sim_d = ss(A_SIM_D, B_SIM_D, C_SIM_D, D_SIM_D);
    pzmap(sys_sim_d,'r');
    hold on
end
toc

pzmap(sys_real,'b');
hold on;