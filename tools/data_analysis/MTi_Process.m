function y = MTi_Process(u, method)
%将MTi获取的陀螺仪数据转换为角加速度
%   u是MTi数据，为单列向量
%   method是差分方法，可选值包括：'center'：中心差分,'back'：后向差分
%   y是角加速度数据，为单列向量，'center'时，长度比u少2，'back'时，长度比u少1

L = length(u);
switch(method)
    case 'center'
        y = zeros(L-2, 1);
        for i = 2:L-1
            y(i-1) = (u(i + 1) - u(i - 1)) / 2;
        end
    case 'back'
        y = diff(u);
end