function M = ARX_LS_SYSID(InputData, OutputData, Fs, TrainLength)
%利用最小二乘法，采用ARX模型，对系统进行辨识。
%   InputData是输入数据
%   OutputData是输出数据
%   Fs是采样频率，单位Hz
%   TrainLength是利用AIC准则得到系统阶次所需的数据长度

testdata=iddata(InputData, OutputData, 1 / Fs);
testdata.inputname = 'MEMS陀螺解算的角加速度';
testdata.outputname = '角加速度计输出的电压';
ze1 = testdata(1 : length(InputData));
ze1 = dtrend(ze1);
%ze1 = idfilt(ze1, 2, 0.1);

%求最合适的结构
NN = struc(1:10, 1:10, 1:10);
V = arxstruc(testdata(1 : TrainLength),testdata(TrainLength + 1 : length(InputData)), NN);
nn = selstruc(V, 'aic');
M = arx(ze1, nn);

figure,plot(ze1);
figure,compare(ze1,M);
figure,resid(M,ze1);