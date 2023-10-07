function yita_ou=daqi(x,y,z,xt,yt,h)%(x,y,z)镜面中心坐标 (xt,yt,h)镜面中心坐标
dis=((x-xt)^2+(y-yt)^2+(z-h)^2)^0.5;
if dis<=1000
    yita_ou=0.99321-0.0001176*dis+1.97*1e-8*dis^2;
else
    yita_ou=exp(-0.0001106*dis);
end

