clear;
clc;
close all;
warning off;
st=[9 10.30 12 13.30 15]; %当地时间
d=[-59 -29 0 31 61 92 122 153 184 214 245 275]; %春分时间之差
fai=39.4; %纬度
h=80; l=6; w=6; o=1;
sa=zeros(12,5);
sb=zeros(12,5);
DNI=zeros(1,60);
x=readmatrix('附件.xlsx','Range','A2:A1746');
y=readmatrix('附件.xlsx','Range','B2:B1746');
z=readmatrix('附件.xlsx','Range','C2:C1746');
N=length(x);
xt=0; yt=0;
m_zb=[x,y,z];
yita_sb=zeros(60,N);
yita_ou=zeros(60,N);
yita_tr=zeros(60,N);
yita_cos=zeros(60,N);
yita_ref=0.92;
A=36;
theta_e=zeros(12,5);
tachang=zeros(12,5);
xying1=zeros(12,5);
xying2=zeros(12,5);
yying1=zeros(12,5);
yying2=zeros(12,5);
for i=1:12
    for j=1:5
        [sa(i,j) , sb(i,j), DNI(1,5*(i-1)+j)]=sun(st(j),d(i),fai);
        theta_e(i,j)=abs(asin(sa(i,j)));
        tachang(i,j)=h*cot(asin(sa(i,j)));
        sita=abs(asin(sb(i,j))-pi);
        xying1(i,j)=tachang(i,j)*cos(sita+pi/2)-sin(sita)*3.5;
        xying2(i,j)=tachang(i,j)*cos(sita+pi/2)+sin(sita)*3.5;
        yying1(i,j)=tachang(i,j)*sin(sita+pi/2)-cos(sita)*3.5;
        yying2(i,j)=tachang(i,j)*sin(sita+pi/2)+cos(sita)*3.5;
        xv=[xying1(i,j), xying2(i,j),0];
        yv=[yying1(i,j),yying2(i,j),0];
        [vl, vr, vn]=guangxian(x,y,z,h,sa(i,j),sb(i,j));
        for k=1:N
             yita_cos(5*(i-1)+j,k)=abs(vl(1)*vn(k,1)+vl(2)*vn(k,2)+vl(3)*vn(k,3));
            [r_zb,f_zb,r_N,f_N]=if_zd1(m_zb(k,:),m_zb,l,w,theta_e(i,j),vr,vn);
            mb_N=vn(k,:);
            yita_sb(5*(i-1)+j,k)=yita_sb1(m_zb(k,:),r_zb,f_zb,mb_N,r_N,f_N,l,w,vl);
            yita_ou(5*(i-1)+j,k)=daqi(m_zb(k,1),m_zb(k,2),m_zb(k,3),xt,yt,h);
            yita_tr(5*(i-1)+j,k)=jieduan(m_zb(k,1),m_zb(k,2),m_zb(k,3),xt,yt,h,vr);
        end
        [in,~]=inpolygon(x,y,xv,yv);
        yita_sb((5*(i-1)+j),in)=0;
    end
end
yita=yita_tr.*yita_ou.*yita_sb.*yita_cos*yita_ref;
psun1=A*(yita)';
psun2=DNI.*sum(psun1);
yita_cosx=(yita_cos)';
yita_sbx=(yita_sb)';
yita_trx=(yita_tr)';
yitax=(yita)';
yita_cos2=sum(yita_cosx)/N;
yita_sb2=sum(yita_sbx)/N;
yita_tr2=sum(yita_trx)/N;
yita2=sum(yitax)/N;
for i=1:5:60
    yita_cos3(o)=(yita_cos2(i)+yita_cos2(i+1)+yita_cos2(i+2)+yita_cos2(i+3)+yita_cos2(i+4))/5;
    yita_sb3(o)=(yita_sb2(i)+yita_sb2(i+1)+yita_sb2(i+2)+yita_sb2(i+3)+yita_sb2(i+4))/5;
    yita_tr3(o)=(yita_tr2(i)+yita_tr2(i+1)+yita_tr2(i+2)+yita_tr2(i+3)+yita_tr2(i+4))/5;
    psun3(o)=(psun2(i)+psun2(i+1)+psun2(i+2)+psun2(i+3)+psun2(i+4))/5;
    yita3(o)=(yita2(i)+yita2(i+1)+yita2(i+2)+yita2(i+3)+yita2(i+4))/5;
    o=o+1;
end
psun=(sum(psun2)/60)/1000;
yitatot=(sum(yita3)/12);
yita_costot=(sum(yita_cos3)/12);
yita_trtot=(sum(yita_tr3)/12);
yita_sbtot=(sum(yita_sb3)/12);
repsun3=psun3/(N*A);
repsun=psun/(N*A)*1000;
for i=1:12
    fprintf('%.d月21日平均光学效率为%f\n',i,yita3(i));
    fprintf('%d月21日平均阴影遮挡效率为%f\n',i,yita_sb3(i));
    fprintf('%d月21日平均余弦遮挡效率为%f\n',i,yita_cos3(i));
    fprintf('%d月21日平均截断效率为%f\n',i,yita_tr3(i));
    fprintf('%d月21日单位面积镜面平均输出热功率为%f\n',i,repsun3(i));
end
fprintf('年平均光学效率为%f\n',yitatot);
fprintf('年平均余弦效率为%f\n',yita_costot);
fprintf('年平均阴影遮挡效率为%f\n',yita_sbtot);
fprintf('年平均截断效率为%f\n',yita_trtot);
fprintf('年平均输出热功率为%f\n',psun);
fprintf('年平均单位面积镜面平均输出热功率为%f\n',repsun);



