clc,clear,warning off;
X1=table2array(readtable('优化结果3.xlsx'));
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
k=floor((600-sum(nL2))/sum(L2));
n=10*k;m=zeros(n,1);
zb=[];l=[];w=[];
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
    theta=sfi;zb1=zeros(m(i),3);l1=zeros(m(i),1);w1=zeros(m(i),1);
    for p=1:m(i)
        zb1(p,1)=r*cos((p-1)*theta);
        zb1(p,2)=r*sin((p-1)*theta);
        zb1(p,3)=X(cs+1,3);
        l1(p)=X(cs+1,1);
        w1(p)=X(cs+1,2);
    end
    zb=[zb;zb1];l=[l;l1];w=[w;w1];
end
r=sqrt(zb(:,1).^2+(zb(:,2)-X1(31)).^2);
m_zb=zb(r<350,:);w=w(r<350,:);l=l(r<350,:);
n_zb=zb(r>=350,:);
scatter3(m_zb(:,1),m_zb(:,2),m_zb(:,3));
hold on
scatter3(n_zb(:,1),n_zb(:,2),n_zb(:,3));

