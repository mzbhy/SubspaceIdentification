function [y, H] = Butterworth_Filter(u, Fs, fp, fs, Rp, As)
%Butterworth Filter.
%   u是原始数据，为单行（或单列）向量
%   Fs是采样频率，单位Hz
%   fs是阻带截止频率，单位Hz
%   fp是通带截止频率，单位Hz
%   Rp是通带最大衰减
%   As是阻带最小衰减
%   y是滤波后数据，为单行（或单列）向量，长度与u相同
%   H是滤波器参数，用于绘制滤波器频率特性
%   常用低通参数：Fs=5000，fp=1，fs=5，Rp=1dB，As=30dB

wp = fp / Fs;
ws = fs / Fs;
[N, wc] = buttord(wp, ws, Rp, As); %计算率波器的阶数和3dB截止频率
[B, A] = butter(N, wc); %计算滤波器系统函数分子分母多项式
[H, W] = freqz(B, A);

y = filter(B, A, u);
