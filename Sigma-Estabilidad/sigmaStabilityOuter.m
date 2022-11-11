
% function sigmaStabilityOuter(sgma)

set(groot,'defaultAxesTickLabelInterpreter','latex');

figure(1)
clf

figure(3)
clf



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

nv1=(E+VC)/(Co*Lo);
nv2=-(Dii*(IL+ILO))/(2*C*Co*Lo);
nv3=((Di^2+Di*Dii)*(E+VC))/(4*C*Co*L*Lo);

nc1 = (1/(2*L))*(E+VC);
nc2 = (1/(2*L))*(((1/(Co*R))*(E+VC))+((Di/(2*C))*(IL+ILO)));
nc3 = (1/(2*L))*(((Di/(2*C*Co*R))*(IL+ILO))+(((2*C+Co*(Dii^2+Di*Dii))/(2*C*Co*Lo))*(E+VC)));
nc4 = (1/(2*L))*((((Dii*(Di+Dii))/(2*C*Co*Lo*R))*(E+VC))+(((Di)/(2*C*Co*Lo))*(IL+ILO)));

d1=1/(Co*R);
d2=(Co*Lo*Di^2 + 2*Co*L*Dii^2 + 4*C*L)/(4*C*Co*L*Lo);
d3=((Lo*Di^2)+(2*L*Dii^2))/(4*C*Co*L*Lo*R);
d4=(Di^2)/(4*C*Co*L*Lo);

kp1=0.056;
ki1=5.4e3;

% [ 1.9, 9.3e+3]


n1=nv1*kp1;
n2=nv1*ki1+nv2*kp1;
n3=nv2*ki1+nv3*kp1;
n4=nv3*ki1;

VoltageNum = [n1 n2 n3 n4];

cpc1 = d1 + (kp1*nc1);
cpc2 = d2 + (kp1*nc2) + (ki1*nc1);
cpc3 = d3 + (kp1*nc3) + (ki1*nc2);
cpc4 = d4 + (kp1*nc4) + (ki1*nc3);
cpc5 = ki1*nc4;

VoltageDen = [1 cpc1 cpc2 cpc3 cpc4 cpc5];

% sigma_inicial = 0;
% sigma_final = 5000;
% inc = 1000;
% sgma = linspace(5000,0,3);
vcolor = [0, 0.6, .4];

sgma = 0:-500:-3000;

for iter = 1:1:length(sgma)
    w=linspace(0,1e5,100000);
    sigma = sgma(iter);
    vcolor(1)=iter/length(sgma);

    % Definicion de s
s = sigma+1i.*w;

N0 = n1*sigma^3 + n2*sigma^2 + n3*sigma + n4;
D0 = sigma^5 + cpc1*sigma^4 + cpc2*sigma^3 + cpc3*sigma^2 + cpc4*sigma + cpc5;


N = n1*s.^3 + n2*s.^2 + n3*s + n4;
D = s.^5 + cpc1*s.^4 + cpc2*s.^3 + cpc3*s.^2 + cpc4*s + cpc5;

% Recta en w=0
Kp_rect = linspace(-.2,.2,1000);
Ki_rect=-sigma*Kp_rect-sigma*(D0/N0);
hold on

% Barrido en frecuencia
vKp=-real(D./N)-((sigma./w).*imag(D./N));
vKi=(w+((sigma^2)./w)).*imag(D./N);
hold on


cruce = InterX([vKp; vKi], [Kp_rect; Ki_rect]);




 figure(1)
 hold on
 
 if sigma==0;
        plot(vKp,vKi,'r--');
        axis([0 .15 0 1300])
            xlabel('$$K_p$$','FontSize', 24 , 'interpreter', 'latex');
    ylabel('$$K_i$$','FontSize', 24 ,  'interpreter', 'latex');
 else
%          axis([0 .08 0 500])

    plot(vKp,vKi)
    axis([0 .15 0 1300])
    hold on
    plot(Kp_rect, Ki_rect)
    axis([0 .15 0 1300])
    plot(cruce(1,:),cruce(2,:),'ko')


 end
end




%  figure(1)
%     plot(vKp,vKi)
%     axis([-0.05 .08 0 500])
%     hold on
%     plot(Kp_rect, Ki_rect)
%     axis([-0.05 .08 0 500])
%     plot(cruce(1,:),cruce(2,:),'ko')
%     
    
% Kp0=Kp;
% 
% if z==0;
% Ki0=[0 0];
% else
% Ki0=-z*Kp0-z*(D0/N0);
% end
% 
% figure(1)
% hold on

% if z==0;
% plot(Kp,Ki,'k')
% axis([0 .2 0 1000])
% else 
%     [-.03 .01 0 250]
% curve1=plot(Kp,Ki,'b')
% axis([0 0.2 0 1000])
% curve2=plot(Kp0,Ki0,'r')
% 
% mask=Ki0<Ki;
% fx=[Kp(mask),fliplr(Kp(mask))];
% fy=[Ki(mask),fliplr(Ki0(mask))];
% fh=patch(fx,fy,[0.3 0.4 0.8])
% alpha(0.3)
% 
% 
% vcolor=[j,0,.2]
% 
% fh=fill(fx,fy,vcolor)
% end


% figure(2)
% hold on


% if z==0;
% 
% plot3(Kp,Ki,z*ones(size(Kp)),'k')
% 
% else 
% plot3(Kp,Ki,z*ones(size(Kp)),'b')
% plot3(Kp0,Ki0,z*ones(size(Kp)),'b')
% 
% end

% 
while (1)
    figure(1)
    [kp,ki]=ginput(1);
 
cpv1 = cpc1;
cpv2 = cpc2 + (kp*n1);
cpv3 = cpc3 + (kp*n2) + (ki*n1);
cpv4 = cpc4 + (kp*n3) + (ki*n2);
cpv5 = cpc5 + (kp*n4) + (ki*n3);
cpv6 = ki*n4;

P=[1 cpv1 cpv2 cpv3 cpv4 cpv5 cpv6];
vpa(roots(P),2)

Gp=tf(VoltageNum,VoltageDen);
Gc=pid(kp,ki);
sys=series(Gc,Gp);
OL=Gp*Gc;

H=1;

Mc=feedback(sys,H);

figure(3)
step(Mc)
% K=vpa([kp ki],2)
% % figure(4)
% % bodeplot(Mc)
% 
% opts = bodeoptions('cstprefs');
% opts.PhaseVisible = 'off';
% opts.FreqUnits = 'Hz';
% 
% figure(4)
% clf
% bodeplot(OL,opts)
% hold on
% bodeplot(Mc,opts)
%     xlabel('Frequency','FontSize', 10 , 'interpreter', 'latex');
%     ylabel('Magnitude (dB)','FontSize', 10 ,  'interpreter', 'latex'); 
%     title('Bode Diagram','FontSize', 15 ,  'interpreter', 'latex')
end

%% Valores maximos
%  - 580.0 + 3.4e+4i
%  - 580.0 - 3.4e+4i
%   - 40.0 + 6.2e+3i
%   - 40.0 - 6.2e+3i
%  - 1.2e+3 + 470.0i
%  - 1.2e+3 - 470.0i
%  
%  
% K =
%  Kp      Ki
% [ 0.049, 61.0]
