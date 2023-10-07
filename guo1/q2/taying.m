function taying=taying(st,d,fai,h)
d=d';
w=(pi/12)*(st-12);
sc=sin((2*pi*d)/365)*sin((2*pi/365)*23.45); %赤纬角
sa=cos(asin(sc))*cos(fai)*cos(w)+sc*sin(fai); %高度角
sb=(sc-sa.*sin(fai))./cos(asin(sa)).*cos(fai); %方位角
taying=h*cot(asin(sa));

