
close all
clear 
clc 

tspan=[0 0.02];

% For nonlinear evaluation
x0=[0,0,0,0,0]';

% Configuration of the ODEs
opt=odeset('Reltol',1e-6,'Abstol',1e-6);

[t,x]=ode113(@AVGII,tspan,x0);

Lambda4=0.1;

EsTheta4=x(:,5)+(-Lambda4*x(:,4));

figure(1)
plot(t,x(:,4))
figure(2)
plot(t,EsTheta4)

function dx=AVGII(t,x)
%% Datos del convertidor 


R=676; % valor reportado, se asume desconocida
fR=2*pi*40; % frecuencia de conmutacion para carga dinamica


%% Estimaciones

Theta4 = (R*square(fR*t,40));
Theta4(Theta4<0) = 338;

%% Estados 
iL  = x(1);
vC  = x(2);
iLo = x(3);
vCo = x(4);

%% Estimadores
xi  = x(5);

%% Lambdas
Lambda4=0.1;

eta=-Lambda4*vCo;

%% Thetas Estimadas
EsTheta4=xi+eta;


%% Descripcion de Ecuaciones de Estado
dx=zeros(5,1);

dx(1) = ((1/(2*L))*(E*(1+D) - vC*(1-D) )); 
dx(2) = ((1/(2*C))*(iL*(1-D) - iLo*(1+D) )); 
dx(3) = ((1/Lo)*(vC*(1+D) + E*D - vCo  ));
dx(4) = (1/Co)*((iLo - vCo*(1/Theta4)));
dx(5) = Lambda4*((1/Co)*(iLo - vCo*EsTheta4));
end