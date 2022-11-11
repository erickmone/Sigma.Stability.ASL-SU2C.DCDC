%% Lazo de Corriente G=(I/D)
% clear
% clc 

% function sigmaStabilityInner(sgma)
figure(1)
clf

% figure(3)
% clf
% figure(4)
% clf

inx=1;


% figure(2)
% clf

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

% sigma_inicial = 0;
% sigma_final = 5000;
% inc = 500;

% sgma = linspace(0,4000,10);
sgma = 0:-500:-5000;
vcolor = [0, 0.6, .4];


for iter = 1:1:length(sgma)
%     w=linspace(0,1e5,200000);
w=0:0.75:1e5;
    sigma = sgma(iter);
    vcolor(1)=iter/length(sgma);
%     pause(5)
% Definicion de s

s = sigma+1i.*w;

N0 = nc1*sigma^3 + nc2*sigma^2 + nc3*sigma + nc4;
D0 = sigma^4 + d1*sigma^3 + d2*sigma^2 + d3*sigma + d4;

N = nc1*s.^3 + nc2*s.^2+nc3*s+nc4;
D = s.^4 + d1*s.^3 + d2*s.^2 + d3*s + d4;

% Recta en w=0
Kp_rect = linspace(0,2,1000);
Ki_rect=-sigma*Kp_rect-sigma*(D0/N0);
hold on

% Barrido en frecuencia
vKp=-real(D./N)-((sigma./w).*imag(D./N));
vKi=(w+((sigma^2)./w)).*imag(D./N);
hold on

cruce = InterX([vKp; vKi], [Kp_rect; Ki_rect]);
% ERR=0.5e-4;

% ax=length(cruce(1))
% bx=length(cruce(2))
% find((abs(vKp-cruce(1))<ERR))
% find((abs(vKi-cruce(2))<ERR))
% abs(vKi-cruce(2))
% as=find((abs(vKp-cruce(1))<ERR)&(abs(vKi-cruce(2))<ERR))

% cruceKp=find((abs(vKp-cruce(1))<ERR));
% 
% c2=cruceKp(1);
% ERR2=1.5e-4;
% cruceKp2=find(vKp>0);
% c1=cruceKp2(1);
% 
% cruceKp3=find(Kp_rect>0);
% c3=cruceKp3(1);
% 
% XR=[vKp(c1:c2) Kp_rect(c3) vKp(c1)]; YR=[vKi(c1:c2) Ki_rect(c3) vKi(c1)];
% figure(1)
% plot(XR,YR)
% axis([0 .2 0 3000])

%     relleno=fill(XR, YR, vcolor, 'EdgeColor', [0 0 0])
%     set(relleno, 'facealpha', .5)
%     
%     




%      cruce = length(cruce);%cruce (:, 1);
%     for n = 1: length(vKp)
%         if vKp(n) > cruce
%             paro = n;
%             break
%         end
%     end
% mask=Ki_rect<vKi;
% fx=[vKp(mask),fliplr(Kp_rect(mask))];
% fy=[vKi(mask),fliplr(Ki_rect(mask))];
% figure(1)
% hold on
% fh=fill(fx,fy,vcolor, 'EdgeColor', [0 0 0])
% set(fh, 'facealpha', .5)
% axis([0 .2 0 3000])
% 
 figure(1)
 hold on
 
 if sigma==0;
        vKp0=vKp;
        vKi0=vKi;
        plot(vKp,vKi,'r--');
        axis([0 .3 0 8000])
            xlabel('$$K_p$$','FontSize', 24 , 'interpreter', 'latex');
    ylabel('$$K_i$$','FontSize', 24 ,  'interpreter', 'latex'); 
 else 
    
%     figure(1)
%     clf
%     
    plot(vKp,vKi)
    
    
%     hinicial=fill(vKp,vKi,'b')
% hfinal=fill(vKp(3446:5605),vKi(3446:5605),'b')
% set(hfinal,'EdgeColor','none')


%     fill(vKp,vKi,'b')
     axis([0 .3 0 8000])
    hold on
    plot(Kp_rect, Ki_rect,'b')
     axis([0 .3 0 8000])
    plot(cruce(1,:),cruce(2,:),'ko')

% %    Estos son los ejes del paper     axis([0 .32 0 8000])

    
%     delta=sgma(2)*(iter-2);
%     yp=8e3-delta/2;
% 
%     text('Interpreter','latex','String','$$\sigma=$$','Position',[0.125 yp],'FontSize',18);
%     text(0.135,yp,[' ',num2str(sgma(iter))],'FontSize',14);

    
    
%     plot(vKp0,vKi0,'r')
 end
 
 %     relleno=fill(vKp(1:paro), vKi(1:paro), vcolor, 'EdgeColor', [0 0 0])
%     set(relleno, 'facealpha', .5)
    
% 
% if z==0;
% Ki0=[0 0];
% else
% Ki0=-z*Kp0-z*(D0/N0);
% 
% end
% 
% figure(1)
% hold on
% 
% if z==0;
% plot(Kp,Ki,'r--');
% axis([0 .2 0 3000])
% else 
% %     [-.03 .01 0 250]
% plot(Kp,Ki,'b')
% axis([0 .2 0 3000])
% 
% plot(Kp0,Ki0,'r')
% axis([0 .2 0 3000])
% 
% %% Mascara
% mask=Ki0<Ki;
% fx=[Kp(mask),fliplr(Kp(mask))];
% fy=[Ki(mask),fliplr(Ki0(mask))];
% % fh=patch(fx,fy,[0.3 0.4 0.8])
% % alpha(0.3)
% 
% n=length(Kp);
% Kp1=[Kp,Kp]; Ki1=[zeros(1,n),Kp(n:-1:1)];
% fill(Kp1,Ki1,'g')
% alpha(0.3)
% 
% % polyout=union(fx,fy);
% % plot(polyout,'g')
% % vcolor=[j,0,.2]
% % 
% % fh=fill(fx,fy,vcolor)
% end
% 
% 
% figure(2)
% hold on
% 
% 
% if z==0;
% 
% plot3(Kp,Ki,z*ones(size(Kp)),'r--')
% 
% else 
% plot3(Kp,Ki,z*ones(size(Kp)),'b')
% fill3(Kp,Ki,z*ones(size(Kp)),'c')
% plot3(Kp0,Ki0,z*ones(size(Kp)),'r')
% fill3(Kp0,Ki0,z*ones(size(Kp)),'c')
% alpha(0.3)
% end
figure(2)
bs=-sigma.*(D0./N0);
plot(-sigma,bs,'.')
end

while (1)
    figure(1)
    [kp,ki]=ginput(1);
    
cpc1 = d1 + (kp*nc1);
cpc2 = d2 + (kp*nc2) + (ki*nc1);
cpc3 = d3 + (kp*nc3) + (ki*nc2);
cpc4 = d4 + (kp*nc4) + (ki*nc3);
cpc5 = ki*nc4;

P=[1 cpc1 cpc2 cpc3 cpc4 cpc5];
vpa(roots(P),2)
roots(P)


Gp=tf(CurrentNum,GenericDen);
Gc=pid(kp,ki);
sys=series(Gc,Gp);
H=1;
OL=Gp*Gc;
Mc=feedback(sys,H);

figure(3)
step(Mc)
K=vpa([kp ki],2)

opts = bodeoptions('cstprefs');
opts.PhaseVisible = 'off';
opts.FreqUnits = 'Hz';

figure(4)
bodeplot(OL,opts)
hold on
bodeplot(Mc,opts)

figure(4)
bodeplot(Mc)

%% Ganancias Kp, Ki con polos del lazo interno
% [ 0.1, 240.0]
%  - 960.0 + 3.5e+4i
%  - 960.0 - 3.5e+4i
%            -3.0e+4
%            -5.9e+3
%            -1.8e+3

% [ 0.88, 2.2e+3]
% 
%             -3.1e+5
%  - 2.0e+3 + 3.4e+4i
%  - 2.0e+3 - 3.4e+4i
%             -3.9e+3
%             -2.4e+3

end
% 
% 
% %% Valores finales 
% %  - 510.0 + 3.4e+4i
% %  - 510.0 - 3.4e+4i
% %  - 670.0 + 6.3e+3i
% %  - 670.0 - 6.3e+3i
% %            -1.1e+3
% %  
% %  
% % K =
% %  
% % [ 1.5e-3, 33.0]