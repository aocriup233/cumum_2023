clc,clear,warning off;
x0=[2*ones(1,30),0];
lb=[2*ones(1,30),-346.5];
ub=[8*ones(1,10),8*ones(1,10),6*ones(1,10),346.5];
[x,y]=ga(@mb2,31,[],[],[],[],lb,ub,@nlon)

function y=mb2(X1)
X=[X1(1:10)',X1(11:20)',X1(21:30)'];
L1=zeros(10,1);
L2=zeros(10,1);
nL2=zeros(9,1);
cost=sqrt(1-0.2671^2);
tant=0.2671/cost;
A=X(:,1).*X(:,2);
for i=1:10
    L1(i)=2*X(i,1);
end
for i=1:10
    L2(i)=X(i,2)/0.2671;
end
for i=1:9
    nL2(i)=X(i,2)/(2*0.2671)+X(i+1,2)/2*0.2671+((X(i+1,3)-X(i,3))+X(i+1,2)/2*cost)/tant;
end
k=floor((600-sum(nL2))/sum(L2));
n=10*k;m=zeros(n-1,1);
zb=[];nA=[];
for i=1:n-1
    cs=floor(i/k);
    ys=mod(i,k);
    sumL2=0;sumnL2=0;
    if cs>0&&cs<10
       for j=1:cs
           sumL2=sumL2+k*L2(j);
           sumnL2=sumnL2+k*nL2(j);
       end
    end
    sumL2=sumL2+(ys-1)*L2(cs+1);
    m(i)=floor((100+sumL2+sumnL2)*2*pi/L1(cs+1));
    r=100+sumL2+sumnL2;
    sfi=L1(cs+1)/r;
    theta=sfi;zb1=zeros(m(i),3);nA1=zeros(m(i),1);
    for p=1:m(i)
        zb1(p,1)=r*cos((p-1)*theta);
        zb1(p,2)=r*sin((p-1)*theta);
        zb1(p,3)=X(cs+1,3);
        nA1(p)=A(cs+1);
    end
    zb=[zb;zb1];nA=[nA;nA1];
end

st=[9 10.30 12 13.30 15]; %当地时间
d=[-59 -29 0 31 61 92 122 153 184 214 245 275]; %春分时间之差
fai=39.4; %纬度
h=80; l=6; w=6;
sa=zeros(12,5);
sb=zeros(12,5);
DNI=zeros(1,60);
m_zb=[X(:,1),X(:,2),X(:,3)];
r=sqrt(m_zb(:,1).^2+(m_zb(:,2)-X1(31)).^2);
m_zb=m_zb(r<350,:);
nA=nA(r<350);
N=length(m_zb(:,1));
xt=0; yt=0;
x=m_zb(:,1);y=m_zb(:,2);z=m_zb(:,3);
% yita_sb=zeros(60,N);
yita_ou=zeros(60,N);
yita_tr=zeros(60,N);
yita_cos=zeros(60,N);
yita_ref=0.92;
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
psun1=nA.*(yita)';
psun2=DNI.*sum(psun1);
y=-sum(psun2)/60/sum(nA);
end
function [c,ceq]=nlon(X1)
X=[X1(1:10)',X1(11:20)',X1(21:30)'];
L1=zeros(10,1);
L2=zeros(10,1);
nL2=zeros(9,1);
cost=sqrt(1-0.2671^2);
tant=0.2671/cost;
for i=1:10
    L1(i)=2*X(i,1);
end
for i=1:10
    L2(i)=X(i,2)/0.2671;
end
for i=1:9
    nL2(i)=X(i,2)/(2*0.2671)+X(i+1,2)/2*0.2671+((X(i+1,3)-X(i,3))+X(i+1,2)/2*cost)/tant;
end
c=zeros(29,1);
for i=1:10
    c(i)=5+X(i,2)-L1(i);
    c(10+i)=5+X(i,2)-L2(i);
end
for i=1:9
    c(20+i)=5+X(i,2)-nL2(i);
end
ceq=[];
end