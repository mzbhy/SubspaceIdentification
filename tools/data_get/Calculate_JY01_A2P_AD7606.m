clearvars -EXCEPT VarName1 VarName2
%%%%%删除NAN数据%%%%%
NaNlocation = find(isnan(VarName1));
VarName1(NaNlocation)=[];
VarName2(NaNlocation)=[];


%%%%采样信息设置%%%%%
T = 0.0005; % 采样时间
Fs = 1 / T;              % 采样频率
L = length(VarName1);     % 信号长度
time = T * (0:1:L-1);


%%%%%AD原始数据（JY01）处理%%%%%
for i = 1 : L
    if(VarName1(i) > 32768)
        VarName1(i) = VarName1(i)-65536;
    end
end
Vout = VarName1 * 5 / 32768 * 749 / 249;
beta = Vout * 1000 / 0.5;

hold on
plot(time, beta,'r');
title('角加速度信号');

%%%%%AD采用巴特沃斯滤波器%%%%%
Rp = 1;As = 30;
fp = 20;fs = 30;
wp = fp / Fs;
ws = fs / Fs;
[N, wc] = buttord(wp,ws,Rp,As);%计算率波器的阶数和3dB截止频率
[B, A] = butter(N, wc);%计算滤波器系统函数分子分母多项式
[H, W] = freqz(B, A);
%滤波器特性曲线
figure
plot(W, abs(H))
title('滤波器频率特性');
beta_filter = filter(B, A, beta);
figure
plot(time, beta_filter);
title('JY01滤波结果');


%%%%%计算JY01的fft%%%%%
NFFT = 2^nextpow2(L);   % 与L相邻的下一个2的指数
Y = fft(beta, NFFT) / L;
f = Fs / 2 * linspace(0, 1, NFFT / 2 + 1) / 10;
figure
plot(f,2*abs(Y(1:NFFT/2+1)))
title('JY01的fft');

figure
A2Poutput = VarName2 / 32768 * 5;
P=((A2Poutput-0.5)*5-10)*249;
beta = Vout  * 1000 / 0.5;
A2P_beta=(A2Poutput-mean(A2Poutput)) * 1000 / 0.05;
plot(A2P_beta,'DisplayName','A2P_beta');hold all;plot(beta,'DisplayName','beta');hold off;