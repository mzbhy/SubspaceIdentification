function y = JY01_ADValue_Process(u, vmax)
%将AD7606获取的JY01数据转换为角加速度
%   u是AD数据，为单行（或单列）向量
%   vmax是AD量程，为5或10
%   y是角加速度数据，维数与u相同

L = length(u);
for i = 1 : L
    if(u(i) > 32768)
        u(i) = u(i)-65536;
    end
end
Vout = u * vmax / 32768 * 749 / 249; %运放对JY01进行了缩放，使其在AD测量范围内
y = Vout * 1000 / 0.5;  %标度因数选为0.5（JY01手册中，但该参数需要标定）