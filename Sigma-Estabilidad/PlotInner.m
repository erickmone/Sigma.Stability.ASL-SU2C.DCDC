%% Información del sistema 
clc
set(groot,'defaultAxesTickLabelInterpreter','latex');

figure(1)
clf

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

% % Con sigma definido positivo
sigma_inicial = 0; 
sigma_final = 4000;
inc = 500;
sgma = sigma_inicial:inc:sigma_final;
% sgma = linspace(0,3500,2);


% Se define la Omega

% Color para la grafica de las fronteras
vcolor = [0, 0.6, .4];

%% En este ciclo se obtienen todas las graficas de sigma
for iter = 1:1:length(sgma) 
    sigma = sgma(iter);
    vcolor(1)=iter/length(sgma);
    w=linspace(0,1e5,10000);

    %% Se define a s(sigma;omega)
    s = -sigma+1i.*w;
    
    N0 = nc1*sigma^3 + nc2*sigma^2 + nc3*sigma + nc4;
    D0 = sigma^4 + d1*sigma^3 + d2*sigma^2 + d3*sigma + d4;

    N = nc1*s.^3 + nc2*s.^2+nc3*s+nc4;
    D = s.^4 + d1*s.^3 + d2*s.^2 + d3*s + d4;
    
    %% Se define las rectas en w=0
    Kp_rect = linspace(-0.04,2,10000);
    Ki_rect = sigma*Kp_rect+sigma*(D0/N0);
    hold on 
    
    %% Sweep de Omega
    vKp=-real(D./N)+((sigma./w).*imag(D./N));
    vKi=(w+((sigma^2)./w)).*imag(D./N);
    hold on
    
    %% Se define la intersección entre las rectas de w=0 y el sweep de w>0
    cruce  = InterX([vKp; vKi], [Kp_rect; Ki_rect]);
    
     if sigma==0
        vKp0=vKp;
        vKi0=vKi;
%         lim=[0 .32 0 10000]
%         xline(0.32,'r--')
%         yline(16000,'r--')
        hold on
 mask=Ki_rect<vKi;
fx=[vKp(mask),fliplr(Kp_rect(mask))];
fy=[vKi(mask),fliplr(Ki_rect(mask))];
str = '#f29e4c';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
hinicial=fill(fx,fy,color, 'EdgeColor', [0 0 0])
%             set(hinicial,'EdgeColor','none')
set(gcf,'units','points','position',[10,10,500,470])
        axis([0 .32 0 10000])
         box on
%         plot(Kp_rect, Ki_rect,'r--')
%         plot(cruce(1,:),cruce(2,:),'ko')
            xlabel('$$k_p$$','FontSize', 24 , 'interpreter', 'latex');
            ylabel('$$k_i$$','FontSize', 24 ,  'interpreter', 'latex'); 
            legend('Stability Boundary',  'interpreter', 'latex')
                   

     else 
         %% Regiones
    if sigma==4000
        str = '#2c699a';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
            hfinal=fill(vKp(3492:4555),vKi(3492:4555),color, 'EdgeColor', [0 0 0]) % Para sigma=4000
set(hfinal,'EdgeColor','k')
    end
    
    if sigma==3500
        str = '#048ba8';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
            hfinal=fill(vKp(3464:4935),vKi(3464:4935),color, 'EdgeColor', [0 0 0])
set(hfinal,'EdgeColor','k')
    end 
           
    if sigma==3000
        str = '#0db39e';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
            hfinal=fill(vKp(3445:5610),vKi(3445:5610),color, 'EdgeColor', [0 0 0])
set(hfinal,'EdgeColor','k')
    end 
    
    if sigma==2500
        str = '#16db93';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
            hfinal=fill(vKp(3432:7390),vKi(3432:7390),color, 'EdgeColor', [0 0 0])
set(hfinal,'EdgeColor','k')
    end 
         
    if sigma==2000
        lim=[0 .32 0 16000];
        str = '#83e377';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
            hfinal=fill(vKp(3440:10000),vKi(3440:10000),color, 'EdgeColor', [0 0 0])
set(hfinal,'EdgeColor','k')
    end 
    
    if sigma==1500
        lim=[0.02 .32 0 16000];
        str = '#b9e769';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
            hfinal=fillout(vKp(3000:10000),vKi(3000:10000),lim,color, 'EdgeColor', [0 0 0])
%             set(hfinal,'EdgeColor','none')
set(hfinal,'EdgeColor','k')
    end 
    
    if sigma==1000
        lim=[0.02 .32 0 16000];
        str = '#efea5a';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
            hfinal=fillout(vKp(3000:10000),vKi(3000:10000),lim,color, 'EdgeColor', [0 0 0])
%             set(hfinal,'EdgeColor','none')
set(hfinal,'EdgeColor','k')
    end
    
    if sigma==500
        lim=[0.02 .32 0 16000];
        str = '#f1c453';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
            hfinal=fillout(vKp(3000:10000),vKi(3000:10000),lim,color, 'EdgeColor', [0 0 0])
%             set(hfinal,'EdgeColor','none')
set(hfinal,'EdgeColor','k')
    end
            axis([0 .32 0 10000])
            hold on
%             plot(Kp_rect, Ki_rect,'b')
%             plot(cruce(1,:),cruce(2,:),'ko')
        
    
    
    
    
      end
end


while (1)
    figure(1)
    [kp,ki]=ginput(1)

    
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

figure(2)

step(Mc)
hold on
K=vpa([kp ki],2)

opts = bodeoptions('cstprefs');
opts.PhaseVisible = 'off';
opts.FreqUnits = 'Hz';

figure(3)
bodeplot(OL,opts)
hold on
bodeplot(Mc,opts)

figure(3)
bodeplot(Mc)
end