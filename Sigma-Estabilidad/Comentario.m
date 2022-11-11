clear 
clc
f=10;
overSampRate=30;
fs=overSampRate*f;
T=0:1/fs:1;
duty_cycle=70;
u=square(2*pi*f*T,duty_cycle);
u=(u+1)/2;

[t,y]=ode45(@onda_cuadrada1,transpose(T),[0 0]);
plot(t,y(:,1))
hold on
plot(t,y(:,2))

function dy = onda_cuadrada1(t,y)
R=1000;
Rl=100;
L=41.02e-3;
C=22.51e-9;
f=10;
overSampRate=30;
fs=overSampRate*f;
duty_cycle=70;
u=square(2*pi*f*t,duty_cycle);
u=(u+1)/2;

dy=zeros(2,1);
dy(1)=-(Rl/L)*y(1)-(1/L)*y(2)+(1/L)*u;
dy(2)=(1/C)*y(1)-(1/(R*C))*y(2);
end