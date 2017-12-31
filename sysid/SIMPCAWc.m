function [Ar,Br,Cr,Dr,sigSIMPCAWc]=SIMPCAWc(y,u,kfp,kff,order)

%--------------------------------------------------------------------------

[ny N]=size(y);
if N<ny
    y=y';
    [ny N]=size(y);
end
[nu N]=size(u);
if N<nu
    u=u';
    [nu N]=size(u);
end 

%--------------------------------------------------------------------------

for i=1:N-kff-kfp+1
    for j=1:kfp+kff
        Y((j-1)*ny+1:j*ny,i)=y(:,j+i-1);
    end
end
for i=1:N-kff-kfp+1
    for j=1:kfp+kff
        U((j-1)*nu+1:j*nu,i)=u(:,j+i-1);
    end
end
%--------------------------------------------------------------------------

Yp=Y(1:kfp*ny,:);
Yf=Y(kfp*ny+1:(kfp+kff)*ny,:);
Up=U(1:kfp*nu,:);
Uf=U(kfp*nu+1:(kfp+kff)*nu,:);
%--------------------------------------------------------------------------
Zp=[Yp;Up];
Zf=[Yf;Uf];
% Zp=[Up;Uf];
ZFF=sqrt((1/(N-kff-kfp+1)))*Zf*Zp'*sqrtm(inv(Zp*Zp'));

[P S V]=svd(ZFF);

% Pr=P(:,1:ny*kff-order);
Pr=P(:,(ny*kff)+order+1:end);
Pry=Pr(1:ny*kff,:);
Pru=Pr(ny*kff+1:ny*(kfp+kff),:);

Gammaf=null(Pry');

Cr=Gammaf(1:ny,:);
Ar=pinv(Gammaf(1:ny*(kff-1),:))*Gammaf(ny+1:ny*kff,:);

Fai=-(Pry');
Psai=Pru';

FFai=zeros((ny*kff-order)*kff,(ny*kff));
for i=1:kff
  FFai((ny*kff-order)*(i-1)+1:(ny*kff-order)*i,1:(kff-i+1)*ny)=Fai(:,(i-1)*ny+1:kff*ny);
  PPsai((i-1)*(ny*kff-order)+1:(ny*kff-order)*i,:)=Psai(:,(i-1)*nu+1:i*nu);
end

Hf1=FFai\PPsai;
DB=pinv([eye(ny)  zeros(ny,order); zeros(ny*(kff-1),ny) Gammaf(1:ny*(kff-1),:)])*Hf1;
Dr=DB(1:ny,:);
Br=DB(ny+1:ny+order,:);


for i=1:kff*(ny+nu)

sigSIMPCAWc(i) = S(i,i);
end

       