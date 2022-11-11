
%% ASL-SU2C I&I

close all
clear 
clc 

tspan=[0 0.02];

% For nonlinear evaluation
x0=[0,0,0,0]';

% Configuration of the ODEs
opt=odeset('Reltol',1e-6,'Abstol',1e-6);

[t,x]=ode45(@AVG,tspan,x0);

% %% Lambdas
% Lambda1=0.1;
% Lambda2=0.5;
% Lambda3=0.3;
% Lambda4=0.2;
% 
% %% Thetas Estimadas
% EsTheta1=x(:,5)+(-Lambda1*x(:,1));
% EsTheta2=x(:,6)+(-Lambda2*x(:,2));
% EsTheta3=x(:,7)+(-Lambda3*x(:,3));
% EsTheta4=x(:,8)+(-Lambda4*x(:,4));

%% Datos del convertidor 
L=223e-6;
Lo=2.34e-3;
C=1e-6;
Co=1e-6;
E=20;
D=0.75;

R=676; % valor reportado, se asume desconocida
fR=2*pi*40; % frecuencia de conmutacion para carga dinamica


%% Valores reportados de parasitos
rL  = 0.046;
rC  = 0.010;
rLo = 0.412;

%% Estimaciones
Theta1 = rL;
Theta2 = rC;
Theta3 = rLo;
Theta4 = R*square(fR*t,40);
Theta4(Theta4<0) = 338;

figure(1)
p1=plot(t,Theta4,'k');
figure(2)
p2=plot(t,x(:,1));
figure(3)
p3=plot(t,x(:,2));
figure(4)
p4=plot(t,x(:,3));
figure(5)
p5=plot(t,x(:,4));

function dx=AVG(t,x)
%% Datos del convertidor 
L=223e-6;
Lo=2.34e-3;
C=1e-6;
Co=1e-6;
E=20;
D=0.75;

R=338; % valor reportado, se asume desconocida
fR=2*pi*40; % frecuencia de conmutacion para carga dinamica


%% Valores reportados de parasitos
rL  = 0.046;
rC  = 0.010;
rLo = 0.412;

%% Estimaciones
Theta1 = rL;
Theta2 = rC;
Theta3 = rLo;

Theta4 = R*square(fR*t,40);
Theta4(Theta4<0) = 338;

%% Estados 
iL  = x(1);
vC  = x(2);
iLo = x(3);
vCo = x(4);

%% Estimadores
% z1  = x(5);
% z2  = x(6);
% z3  = x(7);
% z4  = x(8);

%% Lambdas
% Lambda1=0.1;
% Lambda2=0.5;
% Lambda3=0.3;
% Lambda4=0.2;

%% Thetas Estimadas
% EsTheta1=z1+(-Lambda1*iL);
% EsTheta2=z2+(-Lambda2*vC);
% EsTheta3=z3+(-Lambda3*iLo);
% EsTheta4=z4+(-Lambda4*vCo);


%% Descripcion de Ecuaciones de Estado
dx=zeros(4,1);

dx(1) = ((1/(2*L))*(E*(1+D) - vC*(1-D)   - Theta1*iL)); 
dx(2) = ((1/(2*C))*(iL*(1-D) - iLo*(1+D) - Theta2*vC)); 
dx(3) = ((1/Lo)*(vC*(1+D) + E*D - vCo    - Theta3*iLo));
dx(4) = (1/Co)*((iLo - vCo/R));
% 
% dx(5) = Lambda1*(1/(2*L))*(E*(1+D) - vC*(1-D) - EsTheta1*iL);
% dx(6) = Lambda2*(1/(2*C))*(iL*(1-D) - iLo*(1+D) - EsTheta2*vC);
% dx(7) = Lambda3*(1/Lo)*(vC*(1+D) + E*D - vCo - EsTheta3*iLo);
% dx(8) = Lambda4*(1/Co)*(iLo - EsTheta4*vCo);
end