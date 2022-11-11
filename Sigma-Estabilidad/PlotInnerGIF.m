%% Información del sistema 
clc

h=figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';

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

% Con sigma definido positivo
sigma_inicial = 0; 
sigma_final = 4000;
inc = 500;
sgma = sigma_inicial:inc:sigma_final;

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
    
    % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256);
      
     if sigma==0
        vKp0=vKp;
        vKi0=vKi;
        hold on
        
        mask=Ki_rect<vKi;
        fx=[vKp(mask),fliplr(Kp_rect(mask))];
        fy=[vKi(mask),fliplr(Ki_rect(mask))];
        
        hinicial=fill(fx,fy,[0, 1, .4], 'EdgeColor', [0 0 0])
        set(hinicial, 'facealpha', .1)
        axis([0 .32 0 10000])
        xlabel('$$K_p$$','FontSize', 24 , 'interpreter', 'latex');
        ylabel('$$K_i$$','FontSize', 24 ,  'interpreter', 'latex'); 
        print -painters
        hold on
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
        end
                
  
         
%% Regiones de Estabilidad para Sigma > 0
    if sigma==4000
            hfinal=fill(vKp(3492:4555),vKi(3492:4555),[.5, .2, .4], 'EdgeColor', [0 0 0]) % Para sigma=4000
            set(hfinal, 'facealpha', .4)
            print -painters
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    
    if sigma==3500
            hfinal=fill(vKp(3464:4935),vKi(3464:4935),[0, .7, .1], 'EdgeColor', [0 0 0])
            set(hfinal, 'facealpha', .4)
            print -painters
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
           
    if sigma==3000
            hfinal=fill(vKp(3445:5610),vKi(3445:5610),[.2, .1, .9], 'EdgeColor', [0 0 0])
            set(hfinal, 'facealpha', .4)
            print -painters
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    
    if sigma==2500
            hfinal=fill(vKp(3432:7390),vKi(3432:7390),[.9, 1, .8], 'EdgeColor', [0 0 0])
            set(hfinal, 'facealpha', .4)
            print -painters
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
    end 
         
    if sigma==2000
            lim=[0 .32 0 16000];
            hfinal=fill(vKp(3440:10000),vKi(3440:10000),[0.2, .1, .4], 'EdgeColor', [0 0 0])
            set(hfinal, 'facealpha', .4)
            print -painters
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
    end 
    
    if sigma==1500
            lim=[0.02 .32 0 16000];
            hfinal=fillout(vKp(3000:10000),vKi(3000:10000),lim,[0.9, 0, .4], 'EdgeColor', [0 0 0])
            set(hfinal,'EdgeColor','k', 'facealpha', .4)
            print -painters
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    
    if sigma==1000
            lim=[0.02 .32 0 16000];
            hfinal=fillout(vKp(3000:10000),vKi(3000:10000),lim,[1, 0, .2], 'EdgeColor', [0 0 0])
            set(hfinal,'EdgeColor','k', 'facealpha', .4)
            print -painters
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    
    if sigma==500
            lim=[0.02 .32 0 16000];
            hfinal=fillout(vKp(3000:10000),vKi(3000:10000),lim,[0, 1, .4], 'EdgeColor', [0 0 0])
            set(hfinal,'EdgeColor','k', 'facealpha', .4)
            print -painters
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else
            imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
            axis([0 .32 0 10000])
            hold on
        
    
      end
