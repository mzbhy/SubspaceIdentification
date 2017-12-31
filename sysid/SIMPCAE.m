function [A,B,C,D,sigSIMPCAE]=SIMPCAE(u,y,f,p,n)

%*******************************************

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

%********************************************

for i=1:N-f-p+1
    for j=1:p+f
        Y((j-1)*ny+1:j*ny,i)=y(:,j+i-1);
    end
end

for i=1:N-f-p+1
    for j=1:p+f
        U((j-1)*nu+1:j*nu,i)=u(:,j+i-1);
    end
end

%********************************************

Yp=Y(1:p*ny,:);
Yf=Y(p*ny+1:(p+f)*ny,:);
Up=U(1:p*nu,:);
Uf=U(p*nu+1:(p+f)*nu,:);

%********************************************

ZQR = [Up;Yp;Uf;Yf];

[Q,L]= qr(ZQR',0);
Q=Q';
L=L';
L22 = L((nu+ny)*p+nu*f+1:end,(nu+ny)*p+nu*f+1:end);
Q2 = Q((nu+ny)*p+nu*f+1:end,:);
Ef=L22*Q2;

Yfb=Yf-Ef;
%********************************************

[P S V]=svd([Yfb;Uf]);

%********************************************

Pr=P(:,(ny*f)+n+1:end);

%********************************************

Pry=Pr(1:ny*f,:);
Pru=Pr(ny*f+1:ny*(f+p),:);
% Pre=Pr(ny*(f+p)+1:end,:);

%********************************************

Gammaf=null(Pry');

%********************************************

C=Gammaf(1:ny,:);
A=pinv(Gammaf(1:ny*(f-1),:))*Gammaf(ny+1:ny*f,:);

%********************************************
Fai=-(Pry');
Psai=Pru';

FFai=zeros((ny*f-n)*f,(ny*f));
for i=1:f
  FFai((ny*f-n)*(i-1)+1:(ny*f-n)*i,1:(f-i+1)*ny)=Fai(:,(i-1)*ny+1:f*ny);
  PPsai((i-1)*(ny*f-n)+1:(ny*f-n)*i,:)=Psai(:,(i-1)*nu+1:i*nu);
end

Hf1=FFai\PPsai;
DB=pinv([eye(ny)  zeros(ny,n); zeros(ny*(f-1),ny) Gammaf(1:ny*(f-1),:)])*Hf1;
D=DB(1:ny,:);
B=DB(ny+1:ny+n,:);
%********************************************

for i=1:f*(ny+nu)
sigSIMPCAE(i) = S(i,i);
end

