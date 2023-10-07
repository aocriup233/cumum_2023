function [yita_tr,delta_tot]=jieduan(x,y,z,xt,yt,h,vr) %(x,y,z)镜面中心坐标 (xt,yt,h)镜面中心坐标
dis2=((x-xt)^2+(y-yt)^2+(z-h)^2)^0.5;
delta_sun=2.51*0.001;
delta_s=0.94*0.001;
delta_bq=(2*delta_s)^2;
hr=8;
rr=3.5;
delta_t=0.63*0.001;
ht=hr*(1-vr(3));
ws=2*rr;
delta_ast=((0.5*(ht^2+ws^2))^0.5)/4*dis2*0.001;
delta_tot=((dis2^2)*(delta_sun^2+delta_bq^2+delta_ast^2+delta_t^2))^0.5;
% pd=makedist('Normal','mu',0,'sigma',delta_tot);
% yita_tr=cdf(pd,-3.5);
fun=@(x,y) 1/(2*pi*delta_tot^2)*exp(-(x.^2+y.^2)/(2*delta_tot^2));
polarfun = @(theta,r) fun(r.*cos(theta),r.*sin(theta)).*r;
yita_tr = integral2(polarfun,0,2*pi,0,300);



