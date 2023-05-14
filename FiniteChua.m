clear all;
close all;
clc
dt=0.0001;T=60;
a=-5/7;b=-3/14;beta=15.6;gama=27;
x1=-0.9;x2=-0.15;x3=1.47;
g=0.05;w1=0.5;w2=0.8;bw=0.1;wh=0.5;

Xdy1=0;Xy1=0;Xy2=0;Xz1=0;
Hf1=0;Sf1=0;hthta(:,1)=[0;0;0];


Hf=0;Sf=0;Xy111=0;Xy121=0;Xy210=0;Xy211=0;Xy212=0;Xy220=0;Xy221=0;Xy222=0;
hgam(1)=0;
afa1=0.001;


for i=1:T/dt
    t=i*dt;
    
    fx=a*x1+b*(abs(x1+1)-abs(x1-1));
    thta1=-beta*(1+a);thta2=beta;thta3=-b*beta;
    
    x1=(beta*(-x1+x2-fx))*dt+x1;
    x2=(x1-x2+x3)*dt+x2;
    x3=(-gama*x2)*dt+x3;
    y1=x1+0*sin(t);y2=x2;
    z1=abs(y1+1)-abs(y1-1);
    
    
  x1p(i)=x1;x2p(i)=x2;x3p(i)=x3;
  
  % %F(t,t)   %Kh
Kht=1-exp(-bw*t);dKht=wh-(wh-bw)*exp(-bw*t);
Xdy1=(dKht*y1-wh*Xdy1)*dt+Xdy1;
Xy1=(Kht*y1-wh*Xy1)*dt+Xy1;
Xy2=(Kht*y2-wh*Xy2)*dt+Xy2;
Xz1=(Kht*z1-wh*Xz1)*dt+Xz1;




 
  %K(t,t)
tau=t;
K1t=exp(-w1*(t-tau))*(1-exp(-bw*tau))^2;
dK1t=(w1)-2*(w1-bw)*exp((0-bw)*t)+(w1-2*bw)*exp((0-2*bw)*t);
ddK1t=exp(-w1*t)*(w1^2*exp(w1*tau)-2*(w1-bw)^2*exp((w1-bw)*tau)+(w1-2*bw)^2*exp((w1-2*bw)*tau));
K2t=exp(-w2*(t-tau))*(1-exp(-bw*tau))^2;
dK2t=w2-2*(w2-bw)*exp((-bw)*t)+(w2-2*bw)*exp((-2*bw)*t);
ddK2t=(w2^2-2*(w2-bw)^2*exp((0-bw)*tau)+(w2-2*bw)^2*exp((0-2*bw)*tau));

%[Vky](t)
%Xy110=(K1t*y1-w1*Xy110)*dt+Xy110;
Xy111=(dK1t*y1-w1*Xy111)*dt+Xy111;
%Xy120=(K2t*y1-w2*Xy120)*dt+Xy120;
Xy121=(dK2t*y1-w2*Xy121)*dt+Xy121;

Xy210=(K1t*y2-w1*Xy210)*dt+Xy210;
Xy211=(dK1t*y2-w1*Xy211)*dt+Xy211;
Xy212=(ddK1t*y2-w1*Xy212)*dt+Xy212;

Xy220=(K2t*y2-w2*Xy220)*dt+Xy220;
Xy221=(dK2t*y2-w2*Xy221)*dt+Xy221;
Xy222=(ddK2t*y2-w2*Xy222)*dt+Xy222;

%bianliang
% bS1=y2*(K1t-dK1t)-y1*K1t+Xy111-Xy211-Xy212;
% bS2=y2*(K2t-dK2t)-y1*K2t+Xy121-Xy221-Xy222;
%bV1=-Xy210;bV2=-Xy220;
%wSh=bS1*K2t-bS2*K1t;wVh=bV1*K2t-bV2*K1t;

s2=Xy212+Xy111-Xy211-dK1t*y2;
s3=Xy222+Xy121-Xy221-dK2t*y2;
h2=Xy220-Xy210;
wSh=s2-s3;wVh=h2;
S=wSh*wVh;H=wVh^2;
Hf=(H-g*Hf)*dt+Hf;Sf=(S-g*Sf)*dt+Sf;


%adaptive law
%hthta=[hta1(i);hta2(i);hta3(i)];
Rf=Hf*hgam(i)-Sf;
m(i)=min(abs(eig(Hf)));
% 
if m(i)<10^(-6)
    hgam(i+1)=hgam(i);
 else
     E1=afa1*abs(Rf)^(1/2).*sign(Rf);
       
    u=(inv(Hf))*(-E1);
hgam(i+1)=(inv(Hf)*(S-g*Sf-(H-g*Hf)*hgam(i))+u)*dt+hgam(i);
end
 

%adaptive law
%hthta=[hta1(i);hta2(i);hta3(i)];
Sh=y1*Kht-Xdy1;Vh=[Xy1;Xy2;Xz1];thta=[thta1;thta2;thta3];
S1=Vh*Sh;H1=Vh*Vh';
Hf1=(H1-g*Hf1)*dt+Hf1;Sf1=(S1-g*Sf1)*dt+Sf1;


Rf1=Hf1*hthta(:,i)-Sf1;
m(i)=min(abs(eig(Hf1)));
% 
if m(i)<10^(-12)
    hthta(:,i+1)=hthta(:,i);
 else
     E1=afa1*abs(Rf1).^(1/2).*sign(Rf1);      
    u1=(inv(Hf1))*(-E1);
hthta(:,i+1)=(inv(Hf1)*(S1-g*Sf1-(H1-g*Hf1)*hthta(:,i))+u1)*dt+hthta(:,i);
 end


% 
%h1p(i)=hta1;h2p(i)=hta2;h3p(i)=hta3;
end
for i=1:(T+dt)/dt
ap(i)=a;bp(i)=b;betp(i)=beta;gamp(i)=gama;
end
 % t=dt:dt:T+dt;
t=dt:dt:T+dt;
figure;
plot(t,gamp,'--',t,hgam);
legend('真实值','估计值')
ylabel('$\gamma$','interpreter','latex');
 xlabel('$t$/$s$','interpreter','latex');
axis([0 T,0,28])
  axes('position',[0.55,0.5,0.25,0.25]); 
 plot(t,gamp,'--',t,hgam)
 set(gca,'YLim',[26.5 27.5])
 set(gca,'xLim',[0 T])

figure;
plot(t,betp,'--',t,hthta(2,:))
legend('真实值','估计值')
%set(gca,'YLim',[0 20])
ylabel('$\beta$','interpreter','latex');
 xlabel('$t$/$s$','interpreter','latex');
axis([0 T,-200,inf])
  axes('position',[0.55,0.5,0.25,0.25]); 
 plot(t,betp,'--',t,hthta(2,:))
 set(gca,'YLim',[15.55 15.65])
 set(gca,'xLim',[0 T])



figure;
plot(t,ap,'--',t,-hthta(1,:).*hthta(2,:).^(-1)-1)
legend('真实值','估计值')
%set(gca,'YLim',[-1 -0.5])
ylabel('$a$','interpreter','latex');
 xlabel('$t$/$s$','interpreter','latex');
axis([0 T,-inf,0])
  axes('position',[0.55,0.5,0.25,0.25]); 
 plot(t,ap,'--',t,-hthta(1,:).*hthta(2,:).^(-1)-1)
 set(gca,'YLim',[-0.72 -0.71])
 set(gca,'xLim',[0 T])


figure;
plot(t,bp,'--',t,-hthta(3,:).*hthta(2,:).^(-1))
legend('真实值','估计值')
%set(gca,'YLim',[-0.35 0])
ylabel('$b$','interpreter','latex');
xlabel('$t$/$s$','interpreter','latex');
axis([0 T,-0.3,inf])
  axes('position',[0.55,0.5,0.25,0.25]); 
 plot(t,bp,'--',t,-hthta(3,:).*hthta(2,:).^(-1))
 set(gca,'YLim',[-0.22 -0.21])
 set(gca,'xLim',[0 T])