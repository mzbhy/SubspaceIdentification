function [A, B, C, D] = SIM_Deterministic(y, u, i, n)
%利用子空间方法对输入输出数据进行辨识
%参考Subspace identification for linear system 56页。
%   y是输出数据
%   u是输入数据
%   i是历史数据的数量，也就是分块Hankel矩阵行数的1/2
%   n是系统阶数，可以不指定
%   A,B,C,D是辨识得到的状态空间矩阵

addpath(genpath('../'));

%   判断输入参数，如果n未制定，准备进行阶次辨识
if (nargin < 4)
    n = [];
end

[m, num_u] = size(u);
if(m > num_u)
    u = u';
    [m, num_u] = size(u);
end

[l, num_y] = size(y);
if(l > num_y)
    y = y';
    [l, num_y] = size(y);
end

if(num_u ~= num_y)
    error('输入输出向量长度必须相等');
end

num = num_u;
j = num - 2 * i + 1;

%   计算Hankel矩阵
Y = blkhank(y, 2 * i, j);
U = blkhank(u, 2 * i, j);
Yp = Y(1:i * l, :);
Yf = Y(i * l + 1:2 * i * l, :);
Up = U(1:i * m, :);
Uf = U(i * m + 1:2 * i * m, :);
Wp = [Up; Yp];


%   计算斜投影
R = triu(qr([U; Y]'))'; %计算QR分解，得到下三角矩阵
R = R(1:2 * i * (m + l), 1:2 * i * (m + l)); %将R阵缩减为方阵

mi2 = 2 * m * i;

%   斜投影为Yf/Uf Wp，则
Rp = [R(1:m * i,:);R(2 * m * i + 1:(2 * m + l) * i, :)]; % Wp对应的R阵
Rf = R((2 * m + l) * i + 1:2 * (m + l) * i, :); 	% Yf对应的R阵
Ru = R(m * i + 1:2 * m * i, :); 		% Uf对应的R阵


%   Van Overschee的例程中计算斜投影的部分，当使用该部分时，需要将上面Ru的第二维度改为1;mi2
% % Perpendicular Future outputs 
% Rfp = [Rf(:,1:mi2) - (Rf(:,1:mi2)/Ru)*Ru,Rf(:,mi2+1:2*(m+l)*i)]; 
% % Perpendicular Past
% Rpp = [Rp(:,1:mi2) - (Rp(:,1:mi2)/Ru)*Ru,Rp(:,mi2+1:2*(m+l)*i)]; 
% 
% if (norm(Rpp(:,(2*m+l)*i-2*l:(2*m+l)*i),'fro')) < 1e-10
%     Ob  = (Rfp*pinv(Rpp')')*Rp;
% else
%     Ob = (Rfp/Rpp)*Rp;
% end

%   参考Subspace Identification For Linear Systems 第166页式（6.1）。但是结果与Van
%   Overschee的例程不同，原因未知。
%   此外，书中公式计算斜投影的第三项是C，而例程中是Rc，造成了矩阵维度不同。
Ob1 = Rf * (eye(2 * (m + l) * i) - Ru' * pinv(Ru * Ru') * Ru);
Ob2 = pinv(Rp * (eye(2 * (m + l) * i) - Ru' * pinv(Ru * Ru') * Ru));
Ob = Ob1 * Ob2 * Rp;


%   计算斜投影的SVD分解
[U,S,~] = svd(Ob);
ss = diag(S);

%   阶次辨识
if(isempty(n))
    x = 1:1:l * i;
    figure(gcf);
    subplot;
    bar(x, ss);
    set(gca,'XLim',[0,10.5]);
    %axis([0, length(ss) + 1, 10^(floor(log10(min(ss)))), 10^(ceil(log10(max(ss))))]);
    title('奇异值');
    xlabel('阶次');
    n = input('系统阶次应为？');
end

U1 = U(:, 1:n);
gam  = U1 * diag(sqrt(ss(1:n)));
gamm = U1(1:l * (i - 1), :)*diag(sqrt(ss(1:n)));

gam_per  = U(:,n+1:l*i)'; 		% 求正交补
gamm_inv = pinv(gamm); 			% 求伪逆
A = gamm_inv * gam(l + 1:l * i, :);
C = gam(1:l, :);

%   计算M和L矩阵，为求解B和D做准备
M = gam_per * (R((2 * m + l) * i + 1:2 * (m + l) * i, :) / R(m * i + 1:2 * m * i, :));
L = gam_per;

%   参考Subspace Identification For Linear Systems 第54页
Lhs = zeros(i * (l * i - n), m);
Rhs = zeros(i * (l * i - n), l * i);
for k = 1:i
  Lhs((k - 1) * (l * i - n) + 1:k * (l * i - n), :) = M(:, (k - 1) * m + 1:k * m);
  Rhs((k - 1) * (l * i - n) + 1:k * (l * i - n), 1:(i - k + 1) * l) = L(:, (k - 1) * l + 1:l * i);
end
Rhs = Rhs*[eye(l), zeros(l, n);zeros(l * (i - 1), l), gamm];

%   求解最小二乘
sol = Rhs \ Lhs;

%   得到B和D矩阵
B = sol(l + 1:l + n, :);
D = sol(1 : l, :);
