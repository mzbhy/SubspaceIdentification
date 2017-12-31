function H = blkhank(y,i,j)
%构建分块Hankel矩阵
%   y是数据向量，一般行数为数据的维度，列数为数据集的长度
%   i是Hankel矩阵的行
%   j是Hankel矩阵的列
%   H是构建得到的Hankel矩阵

[l,nd] = size(y);
if nd < l
    y = y';
    [l,nd] = size(y);
end

% 检查维度
if i < 0
    error('blkHank: i should be positive');
end
if j < 0
    error('blkHank: j should be positive');
end
if j > nd-i+1
    error('blkHank: j too big');
end

% 构建分块Hankel矩阵
H=zeros(l*i,j);
for k=1:i
	H((k-1)*l+1:k*l,:)=y(:,k:k+j-1);
end
