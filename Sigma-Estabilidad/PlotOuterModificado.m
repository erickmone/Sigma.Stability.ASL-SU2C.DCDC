
figure(1)
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

sigma_inicial = 0; 
sigma_final = 3000;
inc = 500;
sgma = sigma_inicial:inc:sigma_final;
vcolor = [0, 0.6, .4];


for iter = 1:1:length(sgma)
    w=linspace(-0.001,1e5,10000);
    sigma = sgma(iter);
    vcolor(1)=iter/length(sgma);

    % Definicion de s
s = -sigma+1i.*w;

N0 = n1*sigma^3 + n2*sigma^2 + n3*sigma + n4;
D0 = sigma^5 + cpc1*sigma^4 + cpc2*sigma^3 + cpc3*sigma^2 + cpc4*sigma + cpc5;


N = n1*s.^3 + n2*s.^2 + n3*s + n4;
D = s.^5 + cpc1*s.^4 + cpc2*s.^3 + cpc3*s.^2 + cpc4*s + cpc5;

% Recta en w=0
Kp_rect = linspace(-.4,.09,4000);
Ki_rect=sigma*Kp_rect+sigma*(D0/N0);
hold on

% Barrido en frecuencia
vKp=-real(D./N)+((sigma./w).*imag(D./N));
vKi=(w+((sigma^2)./w)).*imag(D./N);
hold on


cruce = InterX([vKp; vKi], [Kp_rect; Ki_rect]);

% plot(vKp,vKi)
% hold on
% plot(Kp_rect,Ki_rect)
% plot(cruce(1,:),cruce(2,:),'ko')
% axis([0 .12 0 1400])

if sigma==0
        str = '#f29e4c';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        hinicial=fill(vKp(500:4000),vKi(500:4000),color);
        axis([0 .12 0 1400])
        set(gcf,'units','points','position',[10,10,500,470])
        xlabel('$$k_p^{v}$$','FontSize', 24 , 'interpreter', 'latex');
        ylabel('$$k_i^{v}$$','FontSize', 24 ,  'interpreter', 'latex');
        legend('Stability Boundary','Location','best','interpreter', 'latex')
        hold on 
        box on
else 
    if sigma==500
        x=vKp(301:2280);
        y=vKi(301:2280);

        x1=Kp_rect(1:1980);
        y1=Ki_rect(1:1980);  
%         plot(x,y)
%         hold on
%         plot(x1,y1)
        
        mask=y>y1
        fx = [x1(mask), fliplr(x(mask))];
        fy = [y1(mask), fliplr(y(mask))];
         str = '#f1c453';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        fh = fill(fx,fy,color);
    end
    
    if sigma==1000
        x2=vKp(335:2286);
        y2=vKi(335:2286);

        x3=Kp_rect(1:1952);
        y3=Ki_rect(1:1952);  
%         plot(x2,y2)
%         hold on
%         plot(x3,y3)
        
        mask=y2>y3
        fx = [x3(mask), fliplr(x2(mask))];
        fy = [y3(mask), fliplr(y2(mask))];
         str = '#efea5a';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        fh = fill(fx,fy,color);
    end
    
        if sigma==1500
        x4=vKp(250:2284); 
        y4=vKi(250:2284);

        x5=Kp_rect(1816:3850);
        y5=Ki_rect(1816:3850);  
%         plot(x4,y4)
%         hold on
%         plot(x5,y5)
        
        mask=y4>y5
        fx = [x5(mask), fliplr(x4(mask))];
        fy = [y5(mask), fliplr(y4(mask))];
        % Convert color code to 1-by-3 RGB array (0~1 each)
        str = '#b9e769';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        fh = fill(fx,fy,color);
        end
        
     if sigma==2000
        x6=vKp(348:2268); 
        y6=vKi(348:2268);

        x7=Kp_rect(1917:3837);
        y7=Ki_rect(1917:3837);  
%         plot(x6,y6)
%         hold on
%         plot(x7,y7)
        
        mask=y6>y7
        fx = [x7(mask), fliplr(x6(mask))];
        fy = [y7(mask), fliplr(y6(mask))];
        str = '#83e377';
        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
        fh = fill(fx(1917:3837),fy(1917:3837),color);
        

        end
        
end 





%  if sigma==0
%         hinicial=fill(vKp(500:4000),vKi(500:4000),[1 0.5 1]);
%         axis([0 .12 0 1400])
%         set(gcf,'units','points','position',[10,10,500,500])
%         xlabel('$$k_p$$','FontSize', 24 , 'interpreter', 'latex');
%         ylabel('$$k_i$$','FontSize', 24 ,  'interpreter', 'latex');
%         legend('Stability Boundary','Location','best','interpreter', 'latex')
%         hold on 
%         box on
%  else 
%  if sigma==120
%         x=vKp(276:2278);
%         y=vKi(276:2278);
% 
%         x1=Kp_rect(1:2003);
%         y1=Ki_rect(1:2003);  
% 
%         mask=y>y1
%         fx = [x1(mask), fliplr(x(mask))];
%         fy = [y1(mask), fliplr(y(mask))];
%         fill_color = [1 1 0];
%         fh = fill(fx,fy,fill_color);
% 
%  end
% 
%  if sigma==240
%         x2=vKp(286:2279);
%         y2=vKi(286:2279);
%         x3=Kp_rect(1:1994);
%         y3=Ki_rect(1:1994); 
%         
%         mask=y2>y3
%         fx = [x2(mask), fliplr(x3(mask))];
%         fy = [y2(mask), fliplr(y3(mask))];
%         fill_color = [1 0.6 0];
%         fh = fill(fx,fy,fill_color);
% 
%  end
%  
%  if sigma==360
%         x4=vKp(296:2281);
%         y4=vKi(296:2281);
%         x5=Kp_rect(1:1986);
%         y5=Ki_rect(1:1986); 
%         
%         mask=y4>y5
%         fx = [x4(mask), fliplr(x5(mask))];
%         fy = [y4(mask), fliplr(y5(mask))];
%         fill_color = [1 0.4 0.4];
%         fh = fill(fx,fy,fill_color);
%         
%         
%        
%  end
% end
end