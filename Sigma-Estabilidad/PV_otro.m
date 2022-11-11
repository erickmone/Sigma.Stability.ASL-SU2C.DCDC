clc
clear all

syms V I
%% SRC, also referred as Standard Test Conditions 
%(STC @ 1000W/m^2, 25ºC, AM 1.5)
Pmax11 = 230;
Voc11  = 37.4;
Vmpp11 = 30;
Isc11  = 8.16;
Impp11 = 7.68;

%% Performance at Nominal Operating Condition
% (NOCT @ 800W/m^2, 20ºC, AM 1.5)
Pmax12 = 167;
Voc12  = 33.9;
Vmpp12 = 27.2;
Isc12  = 6.58;
Impp12 = 6.14;

%% Construction
% 60 Cells per Module
% Monocrystaline Silicon Cell Type
% 156x156mm^2 Cell Dimension
Ns=


Isc=2.5; % corriente de corto circuito

Voc_n=0.8391; % voltaje de circuito abierto

Ego=1.103;

alpha_T = 0.0017;

G=[1000,800,600]; %Vector de Radiacion, variando en 1) caso ideal 2) relacion 800/1000 3) relacion 600/1000

R_sh=1000;

R_s=0.02;

T=25+273.15; % Temperatura del medio ambiente

Tr=25+273.15; % Temperatura nominal

If=1.5;   % Factor ideal

Ns=1; % # celdas en serie

G=G';

k=1.38065e-23; % cte de Boltezman

q=1.602e-19; % carga del electron

Vt=If*(k*T/q); 

Rc=2;

n=1;

while n<=length(G)
    
    %% Formulas
    Irs = Isc/((exp((q*Voc_n)/(Ns*k*If*T)))-1);
    
    I_o = Irs*((T/Tr)^3)*exp(((1/Tr)-(1/T))*((q*Ego)/If*k));
    
    I_p = (Isc+alpha_T*(T-Tr))*(G(n,1))/(1000);
        
    f = (I+(-I_p+I_o*((exp((V+I*R_s)/(Vt)))-1)+((V+I*R_s)/(R_sh))));

    v=[0:0.01:1.5]';
    
    [m,~]=size(v);
    
    s1=subs(f,V,v);
    
    w=1;
    
    while w<=m
        s=s1(w,1)==0;
        i_curr=vpasolve(s==0,I);
        i_c(w,n)=double(i_curr);
        w=w+1;
        if i_c(w-1,n)<=0
            break
        end
    end
    vs(1:(w-1),n)=v(1:(w-1),1);
    v_oc(n,1)=vs((w-1),n);
    i_sc(n,1)=i_c(1,n);
    p=vs.*i_c;
    n=n+1;
end

  f2 =(-I+(((Isc+alpha_T*(T-Tr))*(200)/(1000))-I_o*((exp(((I*Rc)+I*R_s)/(Vt)))-1)+(((I*Rc)+I*R_s)/(R_sh))))==0;
  S=solve(f2,I);
i_c(i_c<0) = 0;
p(p<0)=0;
% fprintf('Voltaje de circuito abierto es %10.8f y corriente de corto circuito es %10.8f',v_oc,i_sc)
figure (1)
plot (vs,i_c)
str = sprintf('I(V) Curva Caracteristica para Radiacion diferente');
title(str)
xlabel('V')
ylabel('I')
figure (2)
plot(vs,p)
str = sprintf('P(V) Curva Caracteristica para Radiacion diferente');
title(str)
xlabel('V')
ylabel('P')
