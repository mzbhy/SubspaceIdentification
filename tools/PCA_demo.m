%脚本文件，用于测试PCA

load hald;
[coeff,score,latent,tsquared,explained]= pca(ingredients);
AMean = mean(ingredients);%求数据均值
AStd = std(ingredients);%求数据标准差
B = (ingredients - repmat(AMean,[13 1]))./ repmat(AStd,[13 1]);%求去均值后的原始数据
C = (score(:,1:2))*coeff(:,1:2)' + (score(:,3:4))*coeff(:,3:4)';%测试pca分解:X=PT^{T}+P'T'^{T}