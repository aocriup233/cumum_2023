function [sa, sb, DNI]=sun(st,d,fai) %输入当地时间，时差和纬度
d=d';
w=(pi/12)*(st-12); %太阳时角
g0=1.366; %太阳常数
H=3; %海拔
a= 0.4237 -0.00821*(6-H)^2;
b= 0.5055 + 0.00595*(6.5-H)^2;
c= 0.2711 + 0.01858*(2.5-H)^2;
sc=sin((2*pi*d)/365)*sin((2*pi/365)*23.45); %赤纬角
sa=cos(asin(sc))*cos(fai)*cos(w)+sc*sin(fai); %高度角
sb=(sc-sa.*sin(fai))./cos(asin(sa)).*cos(fai); %方位角
DNI=g0*(a+b*exp(-c./sa)); %单位时间内接收到的太阳辐射能量

