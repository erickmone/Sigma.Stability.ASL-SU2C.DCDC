%% Lazo de Corriente G=(I/D)
% clear
% clc 

function [ven,ved]=modificacion2PIInner(sgma,xm1,xM1)
xm=-10;
xM=10;
figure(1)
clf
ven=[];
ved=[];
inx=1;


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
np=[nc1 nc2 nc3 nc4];

% Denominator
d1=1/(Co*R);
d2=(Co*Lo*Di^2 + 2*Co*L*Dii^2 + 4*C*L)/(4*C*Co*L*Lo);
d3=((Lo*Di^2)+(2*L*Dii^2))/(4*C*Co*L*Lo*R);
d4=(Di^2)/(4*C*Co*L*Lo);
GenericDen=[1 d1 d2 d3 d4];
dp=[1 d1 d2 d3 d4];

% sg=sgma(k);
% sgn=[(-sg)^3; (-sg)^2; -sg; 1];
% sgd=[(-sg)^4; (-sg)^3; (-sg)^2; -sg; 1];



%% Caracterizacion de la Frontera de estabilidad

% sgma = 3400:100:4400;


for iter = 1:1:length(sgma)
    
    w=0.001:0.75:1e5;
    sigma = sgma(iter);
    disp('sgma')
    sigma
    vcolor(1)=iter/length(sgma);
    
    sg=sigma;
    sgn=[-(sg)^3; (sg)^2; -sg; 1];
    sgd=[(sg)^4; -(sg)^3; (sg)^2; -sg; 1];
    ed=dp*sgd;
    en=np*sgn;
    ven(iter)=en;
    ved(iter)=ed;
    bs=sg*(ed/en);
    
    kim=sg*xm+bs;
    kiM=sg*xM+bs;
    
% Definicion de s

s = -sigma+1i.*w;

N0 = nc1*sigma^3 + nc2*sigma^2 + nc3*sigma + nc4;
D0 = sigma^4 + d1*sigma^3 + d2*sigma^2 + d3*sigma + d4;

N = nc1*s.^3 + nc2*s.^2+nc3*s+nc4;
D = s.^4 + d1*s.^3 + d2*s.^2 + d3*s + d4;

% Recta en w=0
Kp_rect = linspace(0,2,1000);
Ki_rect=sigma*Kp_rect+sigma*(D0/N0);
hold on

% Barrido en frecuencia
% vKp=-real(D./N)+((sigma./w).*imag(D./N));
% vKi=(w+((sigma^2)./w)).*imag(D./N);
gs=(s).*(D./N);
vKp=-(1./w).*imag(gs);
vKi=-real(gs)+sigma*vKp;

hold on

cruce = InterX([vKp; vKi], [Kp_rect; Ki_rect]);

 figure(1)
 hold on
 
 if sigma==0;
        plot(vKp,vKi,'r--');
        axis([xm1 xM1 0 8000])
            xlabel('$$K_p$$','FontSize', 24 , 'interpreter', 'latex');
    ylabel('$$K_i$$','FontSize', 24 ,  'interpreter', 'latex'); 
 else 
  
    plot(vKp,vKi)

    hold on
    disp('sigma0')
    sigma
    iter
    plot([xm,xM],[kim,kiM],'m')
     axis([xm1 xM1 0 8000])
    plot(cruce(1,:),cruce(2,:),'ko')

 end

end

% while (1)
%     figure(1)
%     [kp,ki]=ginput(1);
%     
% cpc1 = d1 + (kp*nc1);
% cpc2 = d2 + (kp*nc2) + (ki*nc1);
% cpc3 = d3 + (kp*nc3) + (ki*nc2);
% cpc4 = d4 + (kp*nc4) + (ki*nc3);
% cpc5 = ki*nc4;
% 
% P=[1 cpc1 cpc2 cpc3 cpc4 cpc5];
% roots(P)
% 
% 
% end
