% function [vKp,vKi,kp0,kpF,ki0,kiF,sigmaMax]=sigmaStability(w,sigma)
clc
% close all
figure(1)
clf
figure(2)
clf

%% Caracteristicas del Circuito

% Valores de estado estacionario
IL=70/13;
VC=140;
ILO=10/13;
VCO=260;

% Valores de los componentes del circuito
L=223e-6;
Lo=2.34e-3;
C=1e-6;
Co=1e-6;
R=338;
E=20;
Duty=0.75;
Di=1-Duty;
Dii=1+Duty;

% Numerator
n1=(E+VC)/(Co*Lo);
n2=-(Dii*(IL+ILO))/(2*C*Co*Lo);
n3=((Di^2+Di*Dii)*(E+VC))/(4*C*Co*L*Lo);
VoltageNum=[n1 n2 n3];

% Denominator
d1=1/(Co*R);
d2=(Co*Lo*Di^2 + 2*Co*L*Dii^2 + 4*C*L)/(4*C*Co*L*Lo);
d3=((Lo*Di^2)+(2*L*Dii^2))/(4*C*Co*L*Lo*R);
d4=(Di^2)/(4*C*Co*L*Lo);
VoltageDen=[1 d1 d2 d3 d4];

nc1 = (1/(2*L))*(E+VC);
nc2 = (1/(2*L))*(((1/(Co*R))*(E+VC))+((Di/(2*C))*(IL+ILO)));
nc3 = (1/(2*L))*(((Di/(2*C*Co*R))*(IL+ILO))+(((2*C+Co*(Dii^2+Di*Dii))/(2*C*Co*Lo))*(E+VC)));
nc4 = (1/(2*L))*((((Dii*(Di+Dii))/(2*C*Co*Lo*R))*(E+VC))+(((Di)/(2*C*Co*Lo))*(IL+ILO)));

CurrentNum = [nc1 nc2 nc3 nc4];

%% Caracteristicas de la frontera de estabilidad

sigma_inicial = 0;
sigma_final = 500;
inc = 100;
sgma = sigma_inicial:inc:sigma_final;
vcolor = [0, 0.6, .4];
% endofindex = 0;

for iter = 1:1:length(sgma)
    w=linspace(100,1e4,1000);
    sigma = sgma(iter);
    vcolor(1)=iter/length(sgma);
    
% Definicion de s
s = -sigma+1i.*w;

% Numerador y Denominador
N = n1*s.^2 + n2*s + n3;
D = s.^4 + d1*s.^3 + d2*s.^2 + d3*s + d4;

% Numerador y Denominador
N0 = n1*sigma^2 + n2*sigma + n3;
D0 = sigma^4 + d1*sigma^3 + d2*sigma^2 + d3*sigma + d4;    
    
% Recta en w=0
Kp_rect =  linspace(0,0.001,1000);
Ki_rect=sigma*Kp_rect+sigma*(D0/N0);
hold on
   
% Barrido en frecuencia
vKp=-real(D./N)+((sigma./w).*imag(D./N));
vKi=(w+((sigma^2)./w)).*imag(D./N);
hold on

cruce = InterX([vKp; vKi], [Kp_rect; Ki_rect]);

ERR=1.5e-2;
cruceKp=find((abs(vKp-cruce(1))<ERR)&(abs(vKi-cruce(2))<ERR)&(vKp>0)&(vKp<1e-3));
c2=cruceKp(1)
ERR2=1.5e-4;
cruceKp2=find(vKp>0);
c1=cruceKp2(1);

cruceKp3=find(Kp_rect>0);
c3=cruceKp3(1);

XR=[vKp(c1:c2) Kp_rect(c3) vKp(c1)]; YR=[vKi(c1:c2) Ki_rect(c3) vKi(c1)];

figure(1)
% plot(vKp(c1:c2),vKi(c1:c2))
% axis([0 0.001 0 2])
% hold on
% plot(Kp_rect(c3:end),Ki_rect(c3:end))
% axis([0 0.001 0 2])
plot(XR,YR)
axis([0 0.001 0 1.4])
    relleno=fill(XR, YR, vcolor, 'EdgeColor', [0 0 0])
    set(relleno, 'facealpha', .9)
    xlabel('$$k_p$$','FontSize', 24 , 'interpreter', 'latex');
    ylabel('$$k_i$$','FontSize', 24 ,  'interpreter', 'latex'); 
    legend('Stability Boundary')
    set(gcf,'units','points','position',[10,10,500,470])

%     cruce = cruce (:, 1);
%     for n = 1: length(vKp)
%         if vKp(n) > cruce
%             paro = n;
%             break
%         end
%     end
    
% mask=Ki_rect<vKi;
% fx=[vKp,fliplr(Kp_rect)];
% fy=[vKi,fliplr(Ki_rect)];
% figure(1)
% hold on
% fh=fill(fx,fy,vcolor, 'EdgeColor', [0 0 0])
% set(fh, 'facealpha', .5)
%     axis([0 0.001 0 2])

% hold on
% vColor=[rand rand rand];
% plot(vKp,vKi)
% axis([0 0.001 0 2])
% plot([kp0 kpF],[ki0 kiF],'k','LineWidth',1)
% axis([0 0.001 0 2])

% alpha(0.1)
% fill(vKp,vKi,'r')
% plot(vKp0,vKi0)
% axis([0 0.001 0 2])

% figure (1)
%     variable=fill(vKp(1:endofindex), vKi(1:endofindex), vcolor, 'EdgeColor', [0 0 0]);
%     hold on
%     set(variable, 'facealpha', .8)
%     axis([0 0.001 0 2])
% %     
% figure(1)
% 
% x2=[vKp(1:paro),Kp_rect(1:paro)];
% inBetween=[vKi(1:paro),Ki_rect(1:paro)];
% h=fill(x2, inBetween, vcolor,'EdgeColor', [0 0 0])
% set(h, 'facealpha', .9)

%     figure(1)
%     plot(vKp,vKi)
%     axis([0 0.001 0 2])
%     hold on
%     plot(Kp_rect, Ki_rect,'r')
%     axis([0 0.001 0 2])
%     plot(cruce(1,:),cruce(2,:),'ko')
    
    
%     relleno=fill(vKp(1:paro), vKi(1:paro), vcolor, 'EdgeColor', [0 0 0])
%     set(relleno, 'facealpha', .5)
% Gc=pid(kp,ki);
% sys=series(Gc,Gp);
% H=1;
% 
% Mc=feedback(sys,H);

end

% 
% 
% while (1)
%     figure(1)
%     [kp,ki]=ginput(1);
%     
% 
% Gp=tf(VoltageNum,VoltageDen);
% Gc=pid(kp,ki);
% sys=series(Gc,Gp);
% H=1;
% OL=Gp*Gc;
% 
% Mc=feedback(sys,H);
% opts = bodeoptions('cstprefs');
% opts.PhaseVisible = 'off';
% opts.FreqUnits = 'Hz';
% % figure('Name','Bode Diagram')
% 
% figure(2)
% bodeplot(OL,opts)
% hold on
% bodeplot(Mc,opts)
%     xlabel('Frequency','FontSize', 10 , 'interpreter', 'latex');
%     ylabel('Magnitude (dB)','FontSize', 10 ,  'interpreter', 'latex'); 
%     title('Bode Diagram','FontSize', 15 ,  'interpreter', 'latex')
%     K=vpa([kp ki],2)
% 
% end