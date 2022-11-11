% function [vKp,vKi,kp0,kpF,ki0,kiF,sigmaMax]=sigmaStability(w,sigma)
clc
% close all
figure(1)
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

sigma_inicial = -1;
sigma_final = 500;
inc = 100;
sgma = sigma_inicial:inc:sigma_final;
vcolor = [0, 0.6, .4];
% endofindex = 0;

for iter = 1:1:length(sgma)
    w=linspace(-0.01,1e4,10000);
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
Kp_rect =  linspace(-9e-4,0.001,4000);
Ki_rect=sigma*Kp_rect+sigma*(D0/N0);
hold on
   
% Barrido en frecuencia
vKp=-real(D./N)+((sigma./w).*imag(D./N));
vKi=(w+((sigma^2)./w)).*imag(D./N);
hold on

cruce = InterX([vKp; vKi], [Kp_rect; Ki_rect]);

if sigma==-1
        str = '#f29e4c';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        hinicial=fill(vKp,vKi,color);
        axis([0 0.001 0 1.4])
        set(gcf,'units','points','position',[10,10,500,470])
        xlabel('$$k_p$$','FontSize', 24 , 'interpreter', 'latex');
        ylabel('$$k_i$$','FontSize', 24 ,  'interpreter', 'latex');
        legend('Stability Boundary','Location','best','interpreter', 'latex')
        hold on 
        box on


else
    
    if sigma==99
        x=vKp(3636:7065);
        y=vKi(3636:7065);

        x1=Kp_rect(1:3430);
        y1=Ki_rect(1:3430); 
        
        mask=y>y1
        fx = [x1(mask), fliplr(x(mask))];
        fy = [y1(mask), fliplr(y(mask))];
        
%         plot(x,y)
%         hold on
%         plot(x1,y1)
%         plot(cruce(1,:),cruce(2,:),'ko')
       
        str = '#f1c453';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        fh = fill(fx,fy,color);
    end
    
    if sigma==199
        x2=vKp(3591:6735);
        y2=vKi(3591:6735);

        x3=Kp_rect(1:3145);
        y3=Ki_rect(1:3145); 
        

        
%         plot(x2,y2)
%         hold on
%         plot(x3,y3)
%         plot(cruce(1,:),cruce(2,:),'ko')
        
        mask=y2>y3
        fx = [x3(mask), fliplr(x2(mask))];
        fy = [y3(mask), fliplr(y2(mask))];
         str = '#efea5a';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        fh = fill(fx,fy,color);
    end
    
    if sigma==299
        x4=vKp(3506:6350);
        y4=vKi(3506:6350);
        x5=Kp_rect(1:2845);
        y5=Ki_rect(1:2845); 
        
%         plot(x4,y4)
%         hold on
%         plot(x5,y5)
%         plot(cruce(1,:),cruce(2,:),'ko')
        
        mask=y4>y5
        fx = [x5(mask), fliplr(x4(mask))];
        fy = [y5(mask), fliplr(y4(mask))];
        % Convert color code to 1-by-3 RGB array (0~1 each)
        str = '#b9e769';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        fh = fill(fx,fy,color);
    end
    
    if sigma==399
        x6=vKp(3405:5948);
        y6=vKi(3405:5948);
        x7=Kp_rect(1:2544);
        y7=Ki_rect(1:2544); 
        
%         plot(x6,y6)
%         hold on
%         plot(x7,y7)
%         plot(cruce(1,:),cruce(2,:),'ko')
        
        
        mask=y6>y7
        fx = [x7(mask), fliplr(x6(mask))];
        fy = [y7(mask), fliplr(y6(mask))];
        str = '#83e377';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        fh = fill(fx,fy,color);
    end
    
    if sigma==499
        x8=vKp(3261:5490);
        y8=vKi(3261:5490);
        x9=Kp_rect(1:2230);
        y9=Ki_rect(1:2230); 
        
%         plot(x8,y8)
%         hold on
%         plot(x9,y9)
%         plot(cruce(1,:),cruce(2,:),'ko')
        
        mask=y8>y9
        fx = [x9(mask), fliplr(x8(mask))];
        fy = [y9(mask), fliplr(y8(mask))];
        str = '#16db93';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        fh = fill(fx,fy,color);
    end

% plot(vKp,vKi)
% axis([0 0.001 0 1.4])
% hold on
% plot(Kp_rect,Ki_rect)
% plot(cruce(1,:),cruce(2,:),'ko')

end
end