% Boost Converter TF
clc 

L = 2.7648e-3;
C = 1.666667e-6;
R = 144;
E = 48;
D = 0.6;

Dp = 1-D;
V = E/Dp; % Voltage output desired

Gd3 = 2*V;
Gd2 = (Dp^2)*R;

Gd1 = Gd3/Gd2;
Gd0 = V/Dp;

wy = 2/(R*C);
wz = -((Dp^2)*R)/L;
wo = Dp/(sqrt(L*C));
Q = Dp*R*sqrt(C/L);

numI = [Gd1*((wo^2)/wy) Gd1*(wo^2)];
numV = [Gd0*((wo^2)/wz) Gd0*(wo^2)];
den = [1 wo/Q wo^2];

Gid = tf(numI,den); % minimum fase
Gvd = tf(numV,den); % non-minimum fase

figure(1)
rlocus(Gid,'b',Gvd,'r')
legend({'$G_{id}(s)$','$G_{vd}(s)$'},'Location','best', 'interpreter', 'latex')
hold on

% Coefficients of transfer function Gid

n1 = numI(1);
n0 = numI(2);
d2 = den(1);
d1 = den(2);
d0 = den(3);

% Fronteras de estabilidad

w = linspace(-10000,10000,1000);
s = 1i.*w;

Ng = n1*s + n0;
Dg = d2*s.^2 + d1*s + d0;

% Proposicion Moreno-Mendez-Langarica 2021
vKp = -real(Dg./Ng);
vKi = w.*imag(Dg./Ng);

sigma_inicial = 0; 
sigma_final = 4000;
inc = 1000;
sgma = sigma_inicial:inc:sigma_final;

for iter = 1:1:length(sgma)
    sigma = sgma(iter);
    w = linspace(0,20000,9000);
    sb = -sigma+1i.*w;
    
    Ns0 = n1*sigma + n0;
    Ds0 = d2*sigma^2 + d1*sigma + d0;
    
    Ns01 = n1*sb + n0;
    Ds01 = d2*sb.^2 + d1*sb + d0;
    
    Kp_r = linspace(-0.1,0.1,100);
    Ki_r = sigma*Kp_r + sigma*(Ds0/Ns0);
    
    vKps = -real(Ds01./Ns01)-((sigma./w).*imag(Ds01./Ns01));
    vKis = (w+((sigma^2)./w)).*imag(Ds01./Ns01);
    
    figure(2)
    plot(vKps,vKis)
    xlabel('$$k_p$$','FontSize', 24 , 'interpreter', 'latex');
    ylabel('$$k_i$$','FontSize', 24 ,  'interpreter', 'latex');
    axis([-0.1 0.1 0 1500])
    hold on 
    plot(Kp_r,Ki_r)
      
end

while (1)
    figure(2)
    [kp,ki]=ginput(1)
    
    c3 = d2;
    c2 = (d1 + kp*n1);
    c1 = (d0 + kp*n0 + ki*n1);
    c0 = ki*n0;
    
    e2 = kp*n1;
    e1 = kp*n0 + ki*n1;
    e0 = ki*n0;
    
    Dcl = [c3 c2 c1 c0];
    Ncl = [e2 e1 e0];
    roots(Dcl)
    
    Gcl = tf(Ncl,Dcl)
    figure(3)
    rlocus(Gcl)
    hold on
    
    Gc=pid(kp,ki);
    sys=series(Gc,Gid);
    H=1;
    OL=Gid*Gc;
    Mc=feedback(sys,H);
    
    figure(4)
    bode(Gcl)
    hold on
    
    figure(5)
    step(Gcl)
    hold on
    
end



