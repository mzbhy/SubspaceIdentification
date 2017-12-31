NaNlocation = find(isnan(VarName2));
VarName1(NaNlocation)=[];
VarName2(NaNlocation)=[];
VarName3(NaNlocation)=[];
VarName4(NaNlocation)=[];
VarName5(NaNlocation)=[];

%%%%采样信息设置%%%%%
Fs = 250; % 采样频率
T = 1 / Fs;              % 采样时间
L = length(VarName1);     % 信号长度
time = T * (0:1:L-1);

a0_x=2*9.8*VarName1/32768;
a0_y=2*9.8*VarName2/32768;
a1_x=2*9.8*VarName3/32768;
a2_y=2*9.8*VarName4/32768;

aa=((a2_y-a0_y)-(a1_x-a0_x))*360/2/pi/0.4;
plot(time, aa);