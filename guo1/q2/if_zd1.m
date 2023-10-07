function [r_zb,f_zb,r_N,f_N]=if_zd1(mb_zb,m_zb,l,w,theta_e,Vr,N)
%r_zb:入射遮挡范围坐标,f_zb:反射遮挡范围坐标,m_zb:所有定日镜的坐标,mb_zb:目标定日镜坐标
%l:定日镜长,w:定日镜宽,theta_e:太阳高度角,Vr:主反射光单位向量
L1=2*l;L2=w/sin(theta_e);%判定矩形边长
ux=mb_zb(1)/sqrt(mb_zb(1)^2+mb_zb(2)^2);
uy=mb_zb(2)/sqrt(mb_zb(1)^2+mb_zb(2)^2);
xv=zeros(1,4);yv=zeros(1,4);
xv(2)=mb_zb(1)-L1*uy;
xv(1)=mb_zb(1)+L1*uy;
yv(1)=mb_zb(2)+L1*ux;
yv(2)=mb_zb(2)-L1*ux;
xv(3)=mb_zb(1)+L1*uy+L2*ux;
xv(4)=mb_zb(1)-L1*uy+L2*ux;
yv(3)=mb_zb(2)-L1*ux+L2*uy;
yv(4)=mb_zb(2)+L1*ux+L2*uy;%计算判定矩形坐标
xq=m_zb(:,1);yq=m_zb(:,2);zq=m_zb(:,3);
[in,~]=inpolygon(xq,yq,xv,yv);
r_zb=[xq(in),yq(in),zq(in)];
r_N=N(in,:);
theta_f=asin(abs(Vr(3))/norm(Vr));%计算反射光线与水平面的夹角
L1=2*l;L2=w/sin(theta_f);%判定矩形边长
ux=mb_zb(1)/sqrt(mb_zb(1)^2+mb_zb(2)^2);
uy=mb_zb(2)/sqrt(mb_zb(1)^2+mb_zb(2)^2);
xv=zeros(1,4);yv=zeros(1,4);
xv(1)=mb_zb(1)+L1*uy;
xv(2)=mb_zb(1)-L1*uy;
yv(1)=mb_zb(2)+L1*ux;
yv(2)=mb_zb(2)-L1*ux;
xv(3)=mb_zb(1)+L1*uy+L2*ux;
xv(4)=mb_zb(1)-L1*uy+L2*ux;
yv(3)=mb_zb(2)-L1*ux+L2*uy;
yv(4)=mb_zb(2)+L1*ux+L2*uy;%计算判定矩形坐标
xq=m_zb(:,1);yq=m_zb(:,2);zq=m_zb(:,3);
[in,~]=inpolygon(xq,yq,xv,yv);
f_zb=[xq(in),yq(in),zq(in)];
f_N=N(in,:);







