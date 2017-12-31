function YP = ARX_LS_SYSID_Verify(InputData, OutputData, Fs, M)
%对最小二乘法得到的系统模型进行验证。
%   InputData是输入数据
%   OutputData是输出数据
%   Fs是采样频率，单位Hz
%   M是最小二乘法得到的系统模型
%   YP是InputData在M作用下得到的输出

testdata2=iddata(InputData, OutputData, 1 / Fs);
testdata2.inputname = 'MEMS陀螺解算的角加速度';
testdata2.outputname = '角加速度计输出的电压';
ze2 = testdata2(1 : length(InputData));
ze2 = dtrend(ze2);
%ze2=idfilt(ze2, 2, 0.1);

figure,plot(ze2);
figure,compare(ze2,M);
figure,resid(M,ze2);

%系统仿真
YP = idsim(InputData, M);
