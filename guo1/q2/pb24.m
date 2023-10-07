clc,clear,warning off;
X=[8,5.5645,2];
L2=X(2)/0.2671;
L1=X(1)*2;
n=floor(600/L2);
m=zeros(n,1);
for i=1:n
    m(i)=floor((100+(i-1)*L2)*2*pi/L1);
end
zb=[];
for i=1:n
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
nzb1=zb(sqrt(zb(:,1).^2+zb(:,2).^2)<350,:);
nzb2=zb(sqrt(zb(:,1).^2+zb(:,2).^2)>=350,:);
scatter(nzb1(:,1),nzb1(:,2))
hold on
scatter(nzb2(:,1),nzb2(:,2))