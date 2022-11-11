%% Lazo de Corriente G=(I/D)
% clear
% clc 
close all
clear 
clc

%% Información del sistema 

% Valores en estado estacionario 
IL=70/13;
VC=140;
ILO=10/13;
VCO=260;

% Valores del circuito
L=223e-6;
Lo=2.34e-3;
C=1e-6;
Co=1e-6;
R=338;
E=20;
D=0.75;
Di=1-D;
Dii=1+D;

%% Información de la planta

% Numerator
nc1 = (1/(2*L))*(E+VC);
nc2 = (1/(2*L))*(((1/(Co*R))*(E+VC))+((Di/(2*C))*(IL+ILO)));
nc3 = (1/(2*L))*(((Di/(2*C*Co*R))*(IL+ILO))+(((2*C+Co*(Dii^2+Di*Dii))/(2*C*Co*Lo))*(E+VC)));
nc4 = (1/(2*L))*((((Dii*(Di+Dii))/(2*C*Co*Lo*R))*(E+VC))+(((Di)/(2*C*Co*Lo))*(IL+ILO)));

CurrentNum = [nc1 nc2 nc3 nc4];

% Denominator
d1=1/(Co*R);
d2=(Co*Lo*Di^2 + 2*Co*L*Dii^2 + 4*C*L)/(4*C*Co*L*Lo);
d3=((Lo*Di^2)+(2*L*Dii^2))/(4*C*Co*L*Lo*R);
d4=(Di^2)/(4*C*Co*L*Lo);
GenericDen=[1 d1 d2 d3 d4];

%% Caracterizacion de la Frontera de estabilidad

sgma = 0:-50:-5000; % Sigma para el plot 
% sgma=0:-.05:-20000; % Sigma para ver el limite


for iter = 1:1:length(sgma)

    w=0:0.75:1e5;
    sigma = sgma(iter);

    % Definicion de s

s = sigma+1i.*w;

N0 = nc1*sigma^3 + nc2*sigma^2 + nc3*sigma + nc4;
D0 = sigma^4 + d1*sigma^3 + d2*sigma^2 + d3*sigma + d4;

N = nc1*s.^3 + nc2*s.^2+nc3*s+nc4;
D = s.^4 + d1*s.^3 + d2*s.^2 + d3*s + d4;

% % Recta en w=0
% Kp_rect = linspace(0,2,1000);
% Ki_rect=-sigma*Kp_rect-sigma*(D0/N0);
% hold on
% 
% % Barrido en frecuencia
% vKp=-real(D./N)-((sigma./w).*imag(D./N));
% vKi=(w+((sigma^2)./w)).*imag(D./N);
% hold on
% 
% cruce = InterX([vKp; vKi], [Kp_rect; Ki_rect]);

%  figure(1)
%  hold on
%  
%  if sigma==0
%         vKp0=vKp;
%         vKi0=vKi;
%         plot(vKp,vKi,'r--');
%             axis([0 .3 0 8000])
%             xlabel('$$K_p$$','FontSize', 24 , 'interpreter', 'latex');
%             ylabel('$$K_i$$','FontSize', 24 ,  'interpreter', 'latex'); 
%  else 
% 
%         plot(vKp,vKi)
%      axis([0 .3 0 8000])
%      hold on
%         plot(Kp_rect, Ki_rect,'b')
%      axis([0 .3 0 8000])
%         plot(cruce(1,:),cruce(2,:),'ko')
%  end

 
bs=-sigma.*(D0./N0);
plot(-sigma,bs,'.')

end
