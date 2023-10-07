function [vl, vr, vn]=guangxian(x,y,z,h,sa,sb) %输入镜面中心坐标和吸热器中心离地高度
xl=cos(asin(sa)).*cos(asin(sb)-pi/2);
yl=cos(asin(sa)).*sin(asin(sb)-pi/2);
zl=sa;
vl=[xl ,yl, zl];
vxr=[-x,-y,h-z];
vr=vxr./(x.^2+y.^2+(h-z).^2).^0.5;
m=size(vr);
vxl=repmat(vl,m(1),1);
tv=(vr-vxl);
for i=1:m
    vn(i,1)=tv(i,1)/(tv(i,1)^2+tv(i,2)^2+tv(i,3)^2)^0.5;
    vn(i,2)=tv(i,2)/(tv(i,1)^2+tv(i,2)^2+tv(i,3)^2)^0.5;
    vn(i,3)=tv(i,3)/(tv(i,1)^2+tv(i,2)^2+tv(i,3)^2)^0.5;
end
