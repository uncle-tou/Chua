clear all;
close all;
clc
dt=0.00001;T=5;
x1=-0.9;x2=-0.15;x3=1.47;
gama=27;beta=15.6;a=-5/7;b=-3/14;%选的初始状态及参数
g=0.001;wh=0.1;w=1;%核函数参数 查表
x001=0;x011=0;x021=0;x101=0;x111=0;x121=0;x002=0;x012=0;x022=0;
x102=0;x112=0;x122=0;x202=0;x212=0;x222=0;x003=0;x013=0;x023=0;
s1f=[0;0;0];s2f=0;v1f=[0 0 0; 0 0 0; 0 0 0];v2f=0;
hthe=[0;0;0];hgam=0;
afa1=5;afa2=5;ka1=3/5;ka2=2;%自适应参数 查表





for i=1:T/dt
    t=i*dt;
    f=a*x1+b*(abs(x1+1)-abs(x1-1));
    
%     noise1=0;
%     noise2=0;
    
    noise1=awgn(i-i,10*log10(1000));
    noise2=awgn(i-i,10*log10(1000));
    
   x1=(beta*(-x1+x2-f))*dt+x1;
   x2=(x1-x2+x3)*dt+x2;
   x3=(-gama*x2)*dt+x3;
   y1=noise1+x1;y2=noise2+x2;
   %y1=x1+0.1*rand(1)-0.05;y2=x2+0.1*rand(1)-0.05; 系统建模
   y3=abs(y1+1)-abs(y1-1);
   
   the1=-beta*(1+a);the2=beta;the3=-b*beta;
   

   
%F(t,t)   %Kh
   tau=t;
f0=exp(wh*tau)+exp((wh-2*w)*tau)-2*exp((wh-w)*tau);
df0=wh*exp(wh*tau)+(wh-2*w)*exp((wh-2*w)*tau)-(wh-w)*2*exp((wh-w)*tau);
ddf0=wh^2*exp(wh*tau)+(wh-2*w)^2*exp((wh-2*w)*tau)-(wh-w)^2*2*exp((wh-w)*tau);  

f1=4*exp(wh*tau)-2*exp((wh-w)*tau)-2*exp((wh+w)*tau);
df1=wh*4*exp(wh*tau)-(wh-w)*2*exp((wh-w)*tau)-(wh+w)*2*exp((wh+w)*tau);
ddf1=wh^2*4*exp(wh*tau)-(wh-w)^2*2*exp((wh-w)*tau)-(wh+w)^2*2*exp((wh+w)*tau);

f2=exp(wh*tau)-2*exp((wh+w)*tau)+exp((wh+2*w)*tau);
df2=wh*exp(wh*tau)-(wh+w)*2*exp((wh+w)*tau)+(wh+2*w)*exp((wh+2*w)*tau);
ddf2=wh^2*exp(wh*tau)-(wh+w)^2*2*exp((wh+w)*tau)+(wh+2*w)^2*exp((wh+2*w)*tau);

F00=exp(-wh*t)*f0;F01=exp(-(wh+w)*t)*f1;F02=exp(-(wh+2*w)*t)*f2;
%Kh=exp(-wh*t)*f1+exp(-(wh+w)*t)*f2+exp(-(wh+2*w)*t)*f3;
F10=exp(-wh*t)*df0;F11=exp(-(wh+w)*t)*df1;F12=exp(-(wh+2*w)*t)*df2;
%Kh1=exp(-wh*t)*df1+exp(-(wh+w)*t)*df2+exp(-(wh+2*w)*t)*df3;
F20=exp(-wh*t)*ddf0;F21=exp(-(wh+w)*t)*ddf1;F22=exp(-(wh+2*w)*t)*ddf2;
%Kh2=exp(-wh*t)*ddf1+exp(-(wh+w)*t)*ddf2+exp(-(wh+2*w)*t)*ddf3;


%y1
x001=(F00*y1-(wh+w*0)*x001)*dt+x001;x011=(F01*y1-(wh+w*1)*x011)*dt+x011;x021=(F02*y1-(wh+w*2)*x021)*dt+x021;
V0y1=x001+x011+x021;
x101=(F10*y1-(wh+w*0)*x101)*dt+x101;x111=(F11*y1-(wh+w*1)*x111)*dt+x111;x121=(F12*y1-(wh+w*2)*x121)*dt+x121;
V1y1=x101+x111+x121;
%y2
x002=(F00*y2-(wh+w*0)*x002)*dt+x002;x012=(F01*y2-(wh+w*1)*x012)*dt+x012;x022=(F02*y2-(wh+w*2)*x022)*dt+x022;
V0y2=x002+x012+x022;
x102=(F10*y2-(wh+w*0)*x102)*dt+x102;x112=(F11*y2-(wh+w*1)*x112)*dt+x112;x122=(F12*y2-(wh+w*2)*x122)*dt+x122;
V1y2=x102+x112+x122;
x202=(F20*y2-(wh+w*0)*x202)*dt+x202;x212=(F21*y2-(wh+w*1)*x212)*dt+x212;x222=(F22*y2-(wh+w*2)*x222)*dt+x222;
V2y2=x202+x212+x222;
%y3
x003=(F00*y3-(wh+w*0)*x003)*dt+x003;x013=(F01*y3-(wh+w*1)*x013)*dt+x013;x023=(F02*y3-(wh+w*2)*x023)*dt+x023;
V0y3=x003+x013+x023;




%adaptive law
s1=-V1y1;v1=[V0y1 V0y2 V0y3];s2=V1y1-V1y2+V2y2;v2=-V0y2;
bs1=v1'*s1;bv1=v1'*v1;bs2=v2*s2;bv2=v2^2;

s1f=(-g*s1f+bs1)*dt+s1f;v1f=(-g*v1f+bv1)*dt+v1f;
s2f=(-g*s2f+bs2)*dt+s2f;v2f=(-g*v2f+bv2)*dt+v2f;

R1f=v1f*hthe-s1f;R2f=v2f*hgam-s2f;

m=min(min(abs(eig(v1f))),min(abs(eig(v2f))));
% 
if m<10^(-12)
    hthe=hthe;hgam=hgam;
    t1=t;
else
    E1=-afa1*(abs(R1f).^(ka1)).*sign(R1f)-afa2*(abs(R2f).^(ka2)).*sign(R2f)+(g*v1f-bv1)*hthe-g*s1f+bs1;
    E2=-afa1*abs(R2f)^(ka1)*sign(R2f)-afa2*abs(R2f)^(ka2)*sign(R2f)+(g*v2f-bv2)*hgam-g*s2f+bs2;
    %m1=inv(s2)*(-(-g*s2+Vh)*hta-g*s1+Kh-E1-E2);
   hthe=inv(v1f)*E1*dt+hthe;
   hgam=inv(v2f)*E2*dt+hgam;
end
ap(i)=a;bp(i)=b;betp(i)=beta;gamp(i)=gama;
if abs(hthe(2))>10
ha(i)=-hthe(1)/hthe(2)-1;hb(i)=-hthe(3)/hthe(2);hbet(i)=hthe(2);hg(i)=hgam;
else
    ha(i)=-hthe(1)/10-1;hb(i)=-hthe(3)/10;hbet(i)=hthe(2);hg(i)=hgam;
    
end
%R2fp(i)=R2f;
end
  %t=dt:dt:T+dt;
t=dt:dt:T-dt;

% figure;
% plot(t,-betp+hbet,t,-gamp+hg,t,-ap+ha,t,-bp+hb)%画beta，gama，a，b的图
% xlabel('$t$/$s$','interpreter','latex');
% ylabel('Estimation errors','interpreter','tex');
% axis([0 T,-1,1])


figure;
plot(t,betp,'--',t,hbet)%画实际beta，估计beta
legend('真实值','估计值')
%set(gca,'YLim',[0 20])
ylabel('$\beta$','interpreter','latex');
 xlabel('$t$/$s$','interpreter','latex');
 %灏忓浘
 axes('position',[0.55,0.5,0.25,0.25]); 
 plot(t,betp,'--',t,hbet)
 set(gca,'YLim',[15. 16])
 set(gca,'xLim',[0 T])

figure;
plot(t,gamp,'--',t,hg)%画实际gama，估计gama
legend('真实值','估计值')
%set(gca,'YLim',[0 40])
ylabel('$\gamma$','interpreter','latex');
 xlabel('$t$/$s$','interpreter','latex');
 axes('position',[0.55,0.5,0.25,0.25]); 
 plot(t,gamp,'--',t,hg)
 set(gca,'YLim',[25 30])
 set(gca,'xLim',[0 T])

figure;
plot(t,ap,'--',t,ha)%画实际a，估计a
legend('真实值','估计值')
%set(gca,'YLim',[-1 -0.5])
ylabel('$a$','interpreter','latex');
 xlabel('$t$/$s$','interpreter','latex');
  axes('position',[0.55,0.5,0.25,0.25]); 
 plot(t,ap,'--',t,ha)
 set(gca,'YLim',[-0.75 -0.7])
 set(gca,'xLim',[0 T])


figure;
plot(t,bp,'--',t,hb)%画实际b，估计b
legend('真实值','估计值')
%set(gca,'YLim',[-0.35 0])
ylabel('$b$','interpreter','latex');
 xlabel('$t$/$s$','interpreter','latex');
  axes('position',[0.55,0.5,0.25,0.25]); 
 plot(t,bp,'--',t,hb)
 set(gca,'YLim',[-0.24 -0.2])
 set(gca,'xLim',[0 T])
 





% plot(t,ha1,t,ha2,t,ha3,t,hg)
% set(gca,'YLim',[-0.1 0.1])
% t=dt:dt:T+dt;
% figure;
% subplot(2,2,1); plot(t,z1,t,y);subplot(2,2,2); plot(t,z2,t,y1);
% subplot(2,2,3); plot(t,z3,t,y2);subplot(2,2,4); plot(t,z4,t,y3);
% % figure;
% % plot(t,z1-y,t,z2-y1,t,z3-y2,t,z4-y3)