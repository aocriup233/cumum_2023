function yita_sb=yita_sb1(mb_zb,r_zb,f_zb,mb_N,r_N,f_N,l,w,s_N)
%r_zb:入射遮挡范围坐标,f_zb:反射遮挡范围坐标,mb_zb:目标定日镜坐标,s_N:入射光线方向向量
%l:定日镜长,w:定日镜宽,mb_N:目标定日镜法向量,r_N:入射遮挡定日镜法向量,f_N:反射遮挡定日镜法向量
%蒙特卡洛选点...
    px=-l/2+l*rand(100,1);
    py=-w/2+w*rand(100,1);

[M1,N1]=size(r_zb);[M2,N2]=size(f_zb);
mb_i_zb=zeros(4,3);
sin_beta1=mb_N(1)/sqrt(mb_N(1)^2+mb_N(2)^2);
cos_beta1=-mb_N(2)/sqrt(mb_N(1)^2+mb_N(2)^2);
%计算目标定日镜的三维坐标
mb_i_zb(1,1)=mb_zb(1)-w/2*sqrt(1-mb_N(1)^2)*sin_beta1-l/2*cos_beta1;
mb_i_zb(1,2)=mb_zb(2)+w/2*sqrt(1-mb_N(1)^2)*cos_beta1-l/2*sin_beta1;
mb_i_zb(1,3)=mb_zb(3)+w*mb_N(3)/2;
mb_i_zb(2,1)=mb_zb(1)+w/2*sqrt(1-mb_N(1)^2)*sin_beta1+l/2*cos_beta1;
mb_i_zb(2,2)=mb_zb(2)+w/2*sqrt(1-mb_N(1)^2)*cos_beta1+l/2*sin_beta1;
mb_i_zb(2,3)=mb_zb(3)+w*mb_N(3)/2;
mb_i_zb(3,1)=mb_zb(1)+w/2*sqrt(1-mb_N(1)^2)*sin_beta1+l/2*cos_beta1;
mb_i_zb(3,2)=mb_zb(2)-w/2*sqrt(1-mb_N(1)^2)*cos_beta1+l/2*sin_beta1;
mb_i_zb(3,3)=mb_zb(3)-w*mb_N(3)/2;
mb_i_zb(4,1)=mb_zb(1)-w/2*sqrt(1-mb_N(1)^2)*sin_beta1-l/2*cos_beta1;
mb_i_zb(4,2)=mb_zb(2)-w/2*sqrt(1-mb_N(1)^2)*cos_beta1-l/2*sin_beta1;
mb_i_zb(4,3)=mb_zb(3)-w*mb_N(3)/2;
r_i_zb=zeros(4*M1,N1);f_i_zb=zeros(4*M2,N2);
r_p_zb=zeros(4*M1,N1);f_p_zb=zeros(4*M2,N2);
temp=zeros(4*M1,4);
r_p2_zb=zeros(4*M1,2);f_p2_zb=zeros(4*M2,2);
linel=zeros(4*M1,2);liner=zeros(4*M1,2);
for i=0:M1-1
    %计算每一个入射遮挡定日镜的三维坐标
    sin_betai=r_N(i+1,1)/sqrt(r_N(i+1,1)^2+r_N(i+1,2)^2);
    cos_betai=-r_N(i+1,2)/sqrt(r_N(i+1,1)^2+r_N(i+1,2)^2);
    r_i_zb(4*i+1,1)=r_zb(1)-w/2*sqrt(1-r_N(i+1,1)^2)*sin_betai-l/2*cos_betai;
    r_i_zb(4*i+1,2)=r_zb(2)+w/2*sqrt(1-r_N(i+1,1)^2)*sin_betai-l/2*cos_betai;
    r_i_zb(4*i+1,3)=r_zb(3)+w*r_N(i+1,3)/2;
    r_i_zb(4*i+2,1)=r_zb(1)+w/2*sqrt(1-r_N(i+1,1)^2)*sin_betai+l/2*cos_betai;
    r_i_zb(4*i+2,2)=r_zb(2)+w/2*sqrt(1-r_N(i+1,1)^2)*sin_betai+l/2*cos_betai;
    r_i_zb(4*i+2,3)=r_zb(3)+w*r_N(i+1,3)/2;
    r_i_zb(4*i+3,1)=r_zb(1)+w/2*sqrt(1-r_N(i+1,1)^2)*sin_betai+l/2*cos_betai;
    r_i_zb(4*i+3,2)=r_zb(2)-w/2*sqrt(1-r_N(i+1,1)^2)*sin_betai+l/2*cos_betai;
    r_i_zb(4*i+3,3)=r_zb(3)-w*r_N(i+1,3)/2;
    r_i_zb(4*i+4,1)=r_zb(1)-w/2*sqrt(1-r_N(i+1,1)^2)*sin_betai-l/2*cos_betai;
    r_i_zb(4*i+4,2)=r_zb(2)-w/2*sqrt(1-r_N(i+1,1)^2)*sin_betai-l/2*cos_betai;
    r_i_zb(4*i+4,3)=r_zb(3)-w*r_N(i+1,3)/2;
    %计算每一个入射遮挡定日镜在目标定日镜平面的标准三维投影
    temp(i+1,1)=((mb_i_zb(1,1)-r_i_zb(4*i+1,1)).*mb_N(1)+(mb_i_zb(1,2)-r_i_zb(4*i+1,2)).*mb_N(2)+(mb_i_zb(1,3)-r_i_zb(4*i+1,3)).*mb_N(3))./(s_N*(mb_N)');
    temp(i+1,2)=((mb_i_zb(2,1)-r_i_zb(4*i+2,1)).*mb_N(1)+(mb_i_zb(2,2)-r_i_zb(4*i+2,2)).*mb_N(2)+(mb_i_zb(2,3)-r_i_zb(4*i+2,3)).*mb_N(3))./(s_N*(mb_N)');
    temp(i+1,3)=((mb_i_zb(3,1)-r_i_zb(4*i+3,1)).*mb_N(1)+(mb_i_zb(3,2)-r_i_zb(4*i+3,2)).*mb_N(2)+(mb_i_zb(3,3)-r_i_zb(4*i+3,3)).*mb_N(3))./(s_N*(mb_N)');
    temp(i+1,4)=((mb_i_zb(4,1)-r_i_zb(4*i+4,1)).*mb_N(1)+(mb_i_zb(4,2)-r_i_zb(4*i+4,2)).*mb_N(2)+(mb_i_zb(4,3)-r_i_zb(4*i+4,3)).*mb_N(3))./(s_N*(mb_N)');
    r_p_zb(4*i+1,1)=r_i_zb(4*i+1,1)+(s_N(1).*temp(i+1,1));
    r_p_zb(4*i+1,2)=r_i_zb(4*i+1,1)+(s_N(2).*temp(i+1,1));
    r_p_zb(4*i+1,3)=r_i_zb(4*i+1,1)+(s_N(3).*temp(i+1,1));
    r_p_zb(4*i+2,1)=r_i_zb(4*i+1,1)+(s_N(1).*temp(i+1,2));
    r_p_zb(4*i+2,2)=r_i_zb(4*i+1,1)+(s_N(2).*temp(i+1,2));
    r_p_zb(4*i+2,3)=r_i_zb(4*i+1,1)+(s_N(3).*temp(i+1,2));
    r_p_zb(4*i+3,1)=r_i_zb(4*i+1,1)+(s_N(1).*temp(i+1,3));
    r_p_zb(4*i+3,2)=r_i_zb(4*i+1,1)+(s_N(2).*temp(i+1,3));
    r_p_zb(4*i+3,3)=r_i_zb(4*i+1,1)+(s_N(3).*temp(i+1,3));
    r_p_zb(4*i+4,1)=r_i_zb(4*i+1,1)+(s_N(1).*temp(i+1,4));
    r_p_zb(4*i+4,2)=r_i_zb(4*i+1,1)+(s_N(2).*temp(i+1,4));
    r_p_zb(4*i+4,3)=r_i_zb(4*i+1,1)+(s_N(3).*temp(i+1,4));
    %计算每一个入射遮挡定日镜在目标定日镜平面的平面坐标
    X=[mb_i_zb(2,1)+mb_i_zb(3,1)-2*r_zb(1),mb_i_zb(2,2)+mb_i_zb(3,2)-2*r_zb(2)]/l;
    Y=[mb_i_zb(2,1)+mb_i_zb(1,1)-2*r_zb(1),mb_i_zb(2,2)+mb_i_zb(1,2)-2*r_zb(2)]/l;
    r_p2_zb(4*i+1,1)=(r_p_zb(4*i+1,1)-r_zb(1))*X(1)+(r_p_zb(4*i+1,2)-r_zb(2))*X(2);
    r_p2_zb(4*i+1,2)=(r_p_zb(4*i+1,1)-r_zb(2))*Y(1)+(r_p_zb(4*i+1,2)-r_zb(2))*Y(2);
    r_p2_zb(4*i+2,1)=(r_p_zb(4*i+2,1)-r_zb(1))*X(1)+(r_p_zb(4*i+2,2)-r_zb(2))*X(2);
    r_p2_zb(4*i+2,2)=(r_p_zb(4*i+2,1)-r_zb(2))*Y(1)+(r_p_zb(4*i+2,2)-r_zb(2))*Y(2);
    r_p2_zb(4*i+3,1)=(r_p_zb(4*i+3,1)-r_zb(1))*X(1)+(r_p_zb(4*i+3,2)-r_zb(2))*X(2);
    r_p2_zb(4*i+3,2)=(r_p_zb(4*i+3,1)-r_zb(2))*Y(1)+(r_p_zb(4*i+3,2)-r_zb(2))*Y(2);
    r_p2_zb(4*i+4,1)=(r_p_zb(4*i+1,1)-r_zb(1))*X(1)+(r_p_zb(4*i+4,2)-r_zb(2))*X(2);
    r_p2_zb(4*i+4,2)=(r_p_zb(4*i+1,1)-r_zb(2))*Y(1)+(r_p_zb(4*i+4,2)-r_zb(2))*Y(2);
    %计算每个投影每条边的方程...
    Xl1=[r_p2_zb(4*i+1,1),r_p2_zb(4*i+2,1)];
    Xl2=[r_p2_zb(4*i+1,1),r_p2_zb(4*i+3,1)];
    Xl3=[r_p2_zb(4*i+4,1),r_p2_zb(4*i+2,1)];
    Xl4=[r_p2_zb(4*i+4,1),r_p2_zb(4*i+3,1)];
    Yl1=[r_p2_zb(4*i+1,2),r_p2_zb(4*i+2,2)];
    Yl2=[r_p2_zb(4*i+1,2),r_p2_zb(4*i+3,2)];
    Yl3=[r_p2_zb(4*i+4,2),r_p2_zb(4*i+2,2)];
    Yl4=[r_p2_zb(4*i+4,2),r_p2_zb(4*i+3,2)];
    linel(4*i+1,:)=polyfit(Xl1,Yl1,1);
    linel(4*i+2,:)=polyfit(Xl2,Yl2,1);
    linel(4*i+3,:)=polyfit(Xl3,Yl3,1);
    linel(4*i+4,:)=polyfit(Xl4,Yl4,1);
    %计算与每个投影每条边的交点...
    for j=1:100
        qy1=polyval(linel(4*i+1,:),px(j));
        qy2=polyval(linel(4*i+2,:),px(j));
        qy3=polyval(linel(4*i+3,:),px(j));
        qy4=polyval(linel(4*i+4,:),px(j));
        if (qy1-py(j))*(qy3-py(j))<0 && (qy2-py(j))*(qy4-py(j))<0
            px(j)=10000;
            py(j)=10000;
        end
    end
end
for i=0:M2-1
    %计算每一个反射遮挡定日镜的三维坐标
    sin_betai=f_N(i+1,1)/sqrt(f_N(i+1,1)^2+f_N(i+1,2)^2);
    cos_betai=-f_N(i+1,2)/sqrt(f_N(i+1,1)^2+f_N(i+1,2)^2);
    f_i_zb(4*i+1,1)=f_zb(1)-w/2*sqrt(1-f_N(i+1,1)^2)*sin_betai-l/2*cos_betai;
    f_i_zb(4*i+1,2)=f_zb(2)+w/2*sqrt(1-f_N(i+1,1)^2)*sin_betai-l/2*cos_betai;
    f_i_zb(4*i+1,3)=f_zb(3)+w*f_N(i+1,3)/2;
    f_i_zb(4*i+2,1)=f_zb(1)+w/2*sqrt(1-f_N(i+1,1)^2)*sin_betai+l/2*cos_betai;
    f_i_zb(4*i+2,2)=f_zb(2)+w/2*sqrt(1-f_N(i+1,1)^2)*sin_betai+l/2*cos_betai;
    f_i_zb(4*i+2,3)=f_zb(3)+w*f_N(i+1,3)/2;
    f_i_zb(4*i+3,1)=f_zb(1)+w/2*sqrt(1-f_N(i+1,1)^2)*sin_betai+l/2*cos_betai;
    f_i_zb(4*i+3,2)=f_zb(2)-w/2*sqrt(1-f_N(i+1,1)^2)*sin_betai+l/2*cos_betai;
    f_i_zb(4*i+3,3)=f_zb(3)-w*f_N(i+1,3)/2;
    f_i_zb(4*i+4,1)=f_zb(1)-w/2*sqrt(1-f_N(i+1,1)^2)*sin_betai-l/2*cos_betai;
    f_i_zb(4*i+4,2)=f_zb(2)-w/2*sqrt(1-f_N(i+1,1)^2)*sin_betai-l/2*cos_betai;
    f_i_zb(4*i+4,3)=f_zb(3)-w*f_N(i+1,3)/2;
    %计算每一个反射遮挡定日镜在目标定日镜平面的标准三维投影
    temp(i+1,1)=((mb_i_zb(1,1)-f_i_zb(4*i+1,1))*mb_N(1)+(mb_i_zb(1,2)-f_i_zb(4*i+1,2))*mb_N(2)+(mb_i_zb(1,3)-f_i_zb(4*i+1,3))*mb_N(3))/(s_N*(mb_N)');
    temp(i+1,2)=((mb_i_zb(2,1)-f_i_zb(4*i+2,1))*mb_N(1)+(mb_i_zb(2,2)-f_i_zb(4*i+2,2))*mb_N(2)+(mb_i_zb(2,3)-f_i_zb(4*i+2,3))*mb_N(3))/(s_N*(mb_N)');
    temp(i+1,3)=((mb_i_zb(3,1)-f_i_zb(4*i+3,1))*mb_N(1)+(mb_i_zb(3,2)-f_i_zb(4*i+3,2))*mb_N(2)+(mb_i_zb(3,3)-f_i_zb(4*i+3,3))*mb_N(3))/(s_N*(mb_N)');
    temp(i+1,4)=((mb_i_zb(4,1)-f_i_zb(4*i+4,1))*mb_N(1)+(mb_i_zb(4,2)-f_i_zb(4*i+4,2))*mb_N(2)+(mb_i_zb(4,3)-f_i_zb(4*i+4,3))*mb_N(3))/(s_N*(mb_N)');
    f_p_zb(4*i+1,1)=f_i_zb(4*i+1,1)+(f_N(1)*temp(i+1,1));
    f_p_zb(4*i+1,2)=f_i_zb(4*i+1,1)+(f_N(2)*temp(i+1,1));
    f_p_zb(4*i+1,3)=f_i_zb(4*i+1,1)+(f_N(3)*temp(i+1,1));
    f_p_zb(4*i+2,1)=f_i_zb(4*i+1,1)+(f_N(1)*temp(i+1,2));
    f_p_zb(4*i+2,2)=f_i_zb(4*i+1,1)+(f_N(2)*temp(i+1,2));
    f_p_zb(4*i+2,3)=f_i_zb(4*i+1,1)+(f_N(3)*temp(i+1,2));
    f_p_zb(4*i+3,1)=f_i_zb(4*i+1,1)+(f_N(1)*temp(i+1,3));
    f_p_zb(4*i+3,2)=f_i_zb(4*i+1,1)+(f_N(2)*temp(i+1,3));
    f_p_zb(4*i+3,3)=f_i_zb(4*i+1,1)+(f_N(3)*temp(i+1,3));
    f_p_zb(4*i+4,1)=f_i_zb(4*i+1,1)+(f_N(1)*temp(i+1,4));
    f_p_zb(4*i+4,2)=f_i_zb(4*i+1,1)+(f_N(2)*temp(i+1,4));
    f_p_zb(4*i+4,3)=f_i_zb(4*i+1,1)+(f_N(3)*temp(i+1,4));
    %计算每一个入射遮挡定日镜在目标定日镜平面的平面坐标
    X=[mb_i_zb(2,1)+mb_i_zb(3,1)-2*f_zb(1),mb_i_zb(2,2)+mb_i_zb(3,2)-2*f_zb(2)]/l;
    Y=[mb_i_zb(2,1)+mb_i_zb(1,1)-2*f_zb(1),mb_i_zb(2,2)+mb_i_zb(1,2)-2*f_zb(2)]/l;
    f_p2_zb(4*i+1,1)=(f_p_zb(4*i+1,1)-f_zb(1))*X(1)+(f_p_zb(4*i+1,2)-f_zb(2))*X(2);
    f_p2_zb(4*i+1,2)=(f_p_zb(4*i+1,1)-f_zb(2))*Y(1)+(f_p_zb(4*i+1,2)-f_zb(2))*Y(2);
    f_p2_zb(4*i+2,1)=(f_p_zb(4*i+2,1)-f_zb(1))*X(1)+(f_p_zb(4*i+2,2)-f_zb(2))*X(2);
    f_p2_zb(4*i+2,2)=(f_p_zb(4*i+2,1)-f_zb(2))*Y(1)+(f_p_zb(4*i+2,2)-f_zb(2))*Y(2);
    f_p2_zb(4*i+3,1)=(f_p_zb(4*i+3,1)-f_zb(1))*X(1)+(f_p_zb(4*i+3,2)-f_zb(2))*X(2);
    f_p2_zb(4*i+3,2)=(f_p_zb(4*i+3,1)-f_zb(2))*Y(1)+(f_p_zb(4*i+3,2)-f_zb(2))*Y(2);
    f_p2_zb(4*i+4,1)=(f_p_zb(4*i+1,1)-f_zb(1))*X(1)+(f_p_zb(4*i+4,2)-f_zb(2))*X(2);
    f_p2_zb(4*i+4,2)=(f_p_zb(4*i+1,1)-f_zb(2))*Y(1)+(f_p_zb(4*i+4,2)-f_zb(2))*Y(2);
    %计算每个投影每条边的方程...
    Xr1=[f_p2_zb(4*i+1,1),f_p2_zb(4*i+2,1)];
    Xr2=[f_p2_zb(4*i+1,1),f_p2_zb(4*i+3,1)];
    Xr3=[f_p2_zb(4*i+4,1),f_p2_zb(4*i+2,1)];
    Xr4=[f_p2_zb(4*i+4,1),f_p2_zb(4*i+3,1)];
    Yr1=[f_p2_zb(4*i+1,2),f_p2_zb(4*i+2,2)];
    Yr2=[f_p2_zb(4*i+1,2),f_p2_zb(4*i+3,2)];
    Yr3=[f_p2_zb(4*i+4,2),f_p2_zb(4*i+2,2)];
    Yr4=[f_p2_zb(4*i+4,2),f_p2_zb(4*i+3,2)];
    liner(4*i+1,:)=polyfit(Xr1,Yr1,1);
    liner(4*i+2,:)=polyfit(Xr2,Yr2,1);
    liner(4*i+3,:)=polyfit(Xr3,Yr3,1);
    liner(4*i+4,:)=polyfit(Xr4,Yr4,1);
    for j=1:100
        qy1=polyval(liner(4*i+1,:),px(j));
        qy2=polyval(liner(4*i+2,:),px(j));
        qy3=polyval(liner(4*i+3,:),px(j));
        qy4=polyval(liner(4*i+4,:),px(j));
        if (qy1-py(j))*(qy3-py(j))<0 && (qy2-py(j))*(qy4-py(j))<0
            px(j)=10000;
            py(j)=10000;
        end
    end
end  
    geshu=length(find(px==10000));
    yita_sb=1-(geshu/100);













