function y = FIR_Filter(u, Fs, fp, fs, As)
%Finite Impulse Response Filter，即有限长单位冲激响应滤波器。
%   u是原始数据，为单行（或单列）向量
%   Fs是采样频率，单位Hz
%   fs是阻带截止频率，单位Hz
%   fp是通带截止频率，单位Hz
%   Rp是通带最大衰减
%   As是阻带最小衰减
%   y是滤波后数据，为单行（或单列）向量，长度与u相同
%   常用低通参数：Fs=5000，fp=1，fs=5，Rp=1dB，As=30dB

wp = 2 * fp * pi / Fs; %求归一化数字通带截止频率,
ws = 2 * fs * pi / Fs; %求归一化数字阻带起始频率 
Bt = ws - wp; %求过渡带宽
alpha = 0.5842 * (As - 21)^0.4 + 0.07886 * (As - 21); %计算kaiser窗的控制参数
M = ceil((As - 8) / 2.285 / Bt); %求出滤波器的阶数
wc = (ws + wp) / 2 / pi; %求滤波器的截止频率并关于pi归一化 
hk = fir1(M, wc, kaiser(M + 1, alpha)); %利用 fir1 函数求出滤波器的系数
[Hk, ~] = freqz(hk,1); %计算频率响应
mag = abs(Hk); %求幅频特性
db = 20 * log10(mag / max(mag)); %化为分贝值
db1 = db';
% figure, plot(0 : pi / 511 : pi, db1), grid on
% axis([0, 4.0, -80, 5]), title('fir特性')

y = filter(hk, 1, u);
% figure
% plot(y)
