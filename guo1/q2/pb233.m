clc,clear;warning off;
mypar=parpool;
[x,y]=fmincon(@mb2,[2 2 2 0],[],[],[],[],[2 2 2 -346.5],[8 8 6 346.5],@nlon)

function y=mb2(X)
L2=X(2)/0.2671;
L1=X(1)*2;
n=floor(600/L2);
m=zeros(n,1);
parfor i=1:n
    m(i)=floor((100+(i-1)*L2)*2*pi/L1);
end
zb=[];
parfor i=1:n
    r=100+(i-1)*L2;
    sfi=L1/(100+(i-1)*L2);
    theta=asin(sfi);
    zb1=zeros(m(i),2);
    for j=1:m(i)
        zb1(j,1)=r*cos((j-1)*theta);
        zb1(j,2)=r*sin((j-1)*theta);
    end
    zb=[zb;zb1];
end
z=X(3);
st=[9 10.30 12 13.30 15]; %当地时间
d=[-59 -29 0 31 61 92 122 153 184 214 245 275]; %春分时间之差
fai=39.4; %纬度
h=80; l=6; w=6;
sa=zeros(12,5);
sb=zeros(12,5);
DNI=zeros(1,60);
m_zb=[zb(:,1),zb(:,2)];
r=sqrt(m_zb(:,1).^2+(m_zb(:,2)-X(4)).^2);
m_zb=m_zb(r<350,:);
N=length(m_zb(:,1));
xt=0; yt=0;
m_zb=[m_zb,z*ones(N,1)];
x=m_zb(:,1);y=m_zb(:,2);z=m_zb(:,3);
% yita_sb=zeros(60,N);
yita_ou=zeros(60,N);
yita_tr=zeros(60,N);
yita_cos=zeros(60,N);
yita_ref=0.92;
A=X(1)*X(2);
theta_e=zeros(12,5);
for i=1:12
    for j=1:5
        [sa(i,j), sb(i,j), DNI(1,5*(i-1)+j)]=sun(st(j),d(i),fai);
        theta_e(i,j)=abs(asin(sa(i,j)));
        [vl, vr, vn]=guangxian(x,y,z,h,sa(i,j),sb(i,j));
        yita_cos(5*(i-1)+j,:)=abs(1*vl*(vn)');
        for k=1:N
%             [r_zb,f_zb,r_N,f_N]=if_zd1(m_zb(k,:),m_zb,l,w,theta_e(i,j),vr,vn);
%             mb_N=vn(k,:);
%             yita_sb(5*(i-1)+j,k)=yita_sb1(m_zb(k,:),r_zb,f_zb,mb_N,r_N,f_N,l,w,vl);
            yita_ou(5*(i-1)+j,k)=daqi(m_zb(k,1),m_zb(k,2),m_zb(k,3),xt,yt,h);
            yita_tr(5*(i-1)+j,k)=jieduan(m_zb(k,1),m_zb(k,2),m_zb(k,3),xt,yt,h,vr);
        end
    end
end
yita=yita_tr.*yita_ou.*yita_cos*yita_ref;
psun1=A*(yita)';
psun2=DNI.*sum(psun1);
y=-sum(psun2)/60/A/N+abs(10000*(60000-sum(psun2)/60));
end
function [c,ceq]=nlon(X)
L2=X(2)/0.2671;
L1=X(1)*2;
c(1)=5+X(2)-L1;
c(2)=5+X(2)-L2;
ceq=[];
end