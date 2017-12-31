function[f, ys] = Cal_FFT(u, Fs)
%计算时域序列的FFT并绘图。
%   u是原始数据，为单行（或单列）向量
%   Fs是采样频率，单位Hz

L = length(u);
NFFT = 2^nextpow2(L);   % 与L相邻的下一个2的指数
Y = fft(u, NFFT) / L;
f = Fs / 2 * linspace(0, 1, NFFT / 2 + 1);
ys = 2 * abs(Y(1 : NFFT / 2 + 1));
figure
plot(f, ys)
title('输入信号的单边振幅谱')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

end
  