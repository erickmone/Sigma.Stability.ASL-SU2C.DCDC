%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Junio 2009
%
% Guia 8 - Ejercicio 7
%
%clc
clear all;
warning off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parametros
La=0.02; %[H]
Ra=2;  %[ohm]
kt=0.68; %[Nw m/A]
kw=0.678;  %[V seg/rad]
B=0.02; %[Nw m seg/rad]
J=0.05; %[Nw m seg^2 /rad]
N=0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables de estado:
% w = velocidad angular
% Ia = corriente de armadura
% theta = posición
%
% x=[w Ia theta]'
%
% Matrices del modelo continuo:

Ac=[-Ra/La -kw/La 0;kt/J -B/J 0;0 N 0];
Bc=[1/La 0;0 -N/J; 0 0];
Cc=[1 0 0;0 1 0;0 0 1]; % Asumo que puedo medir todas las variables
Dc=[0 0; 0 0;0 0];

modelo_continuo=ss(Ac,Bc,Cc,Dc);

Uc=[Bc Ac*Bc Ac^2*Bc];
% Si rank(Uc)=3, el sistema es totalmente controlable
rank(Uc) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modelo de estados discreto
%% Periodo de muestreo

Ts=0.003;

% Matrices:

modelo_discreto=c2d(modelo_continuo,Ts);
Ad=modelo_discreto.a;
Bd=modelo_discreto.b;
Cd=modelo_discreto.c;
Dd=modelo_discreto.d;

% Polos del modelo sin realimentar
polos=pole(modelo_discreto); %(z-1)(z-0.7517)(z-0.9843)
transferencia=tf(zpk([],polos,1,Ts));
[num,den] = tfdata(transferencia);
coeficientes=den{1,1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% z^3 + a2 z^2 + a1 z + ao

a0=coeficientes(4)
a1=coeficientes(3)
a2=coeficientes(2)

% Polinomio caracteristico del sistema sin realimentar
% z^3 - 2.736 z^2 + 2.476 z - 0.7399

% Controlabilidad

Ud=[Bd Ad*Bd Ad^2*Bd];
rank(Ud) % El sistema es totalmente controlable

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Realimentacion de estados
%
% u = r - k'x
%
% PO=1%           ts=0.5seg
% PO=100*e^(q)   q=-psi*pi/(1-psi^2)^0.5
% ts=4/(psi*wn)

PO=1;
ts=0.5;
psi=((log(PO/100))^2/(pi^2+(log(PO/100))^2))^0.5
wn=4/(psi*ts)

% s1,2 = -psi*wn ± j wn*(1-psi^2)^0.5 = -8 ± j 5.456

s1=-psi*wn + (wn*(1-psi^2)^0.5)*i;
s2=-psi*wn - (wn*(1-psi^2)^0.5)*i;

% Ubico el tercer polo una decada por encima de la parte real de los polos
% del sistema.
%
% s3 = -80 rad/s
s3 = -80;

% Mapeo al campo discreto z = e^(sT)
% s1->z1  s2->z2 s3->z3
z1=0.976 + 0.01598i;
z2=0.976 - 0.01598i;
z3=0.78;

% Obtención del polinomio caracteristico deseado
% (z-z1)(z-z2)(z-z3)

polos_deseados_discreto=[z1 z2 z3];
transf_deseada_discreta=zpk([],polos_deseados_discreto,1,Ts);
[num_des,den_des] = tfdata(transf_deseada_discreta);
coeficientes_des=den_des{1,1};

% Polinomio caracteristico deseado:
% z^3 - 2.732 z^2 + 2.475 z - 0.7432
% z^3 + a2d   z^2 + a1d   z +  a0d

a0d=coeficientes_des(4) %a0d=-0.7432
a1d=coeficientes_des(3) %a1d= 2.476
a2d=coeficientes_des(2) %a2d=-2.732

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transformacion a un modelo con una sola entrada (Ea)
%  u = [Ea;TL] = f*uf = f*Ea   ==> f=[1;0]

f=[1 ;0];

Af=Ad;
Bf=Bd*f;
Cf=Cd;
Df=Dd;

% Comprobacion de controlabilidad

Uf=[Bf Af*Bf Af^2*Bf];
rank(Uf) % Tiene rango 3 => sigue siendo totalmente controlable

% El vector de realimentacion de estados del modelo con entrada ficticia es
%                  kf' = f^(-1)*k' 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Para calcular kf, utilizo el Modelo Canonico Controlable (MCC). Esto
% simplifica el calculo de los coeficientes del vector de realimentacion
% kcc, el cual queda formado por la diferencia entre los coeficientes del
% polinomio caracteristico sin realimentar y el deseado.
%
% Transformacion al MCC:  Xcc = P X
% P=Ucc*Uf^(-1)           Ucc: matriz controlabilidad del MCC
% kf'=kcc'*P

% Calculo de la matriz de transformacion P
Bcc=[0;0;1];
Acc=[0 1 0; 0 0 1; -a0 -a1 -a2];

Ucc=[Bcc Acc*Bcc Acc^2*Bcc];

P=Ucc*Uf^(-1);

kcc=[a0d-a0; a1d-a1 ; a2d-a2]

kf=((kcc')*P)'

% Vector de realimentacion real:

k=(f*kf')'

k_Ia=k(1,1);
k_omega=k(2,1);
k_theta=k(3,1);

% La salida debe ser 2pi cuando la entrada es 1V:
% k2 es el coeficiente de realimentacion correspondiente a theta:
% Vr - k0*Ia - k1*omega - k2*theta = Ea = 0
%
% Vr = k2*theta ==> 1V = k2 * 2pi
%
% k2=1V/(2pi)
k2=1/(2*pi);

% Ao es la constante utilizada para modificar la ganancia
% k2 = Ao*k_theta
Ao=k2/k_theta;

k_final=k*Ao;

k_Ia=k_final(1,1);
k_omega=k_final(2,1);
k_theta=k_final(3,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Perturbacion TL
% Amplitud maxima de la perturbacion

TLmax=0.5/Ao;

% Calculo un termino proporcional a TL para aplicarlo en forma feed-forward
% y disminuir el efecto sobre la salida.
% En regimen permanente la salida debe ser la misma que sin perturbacion.
%
% Vr - k_theta*theta - k_Ia*Ia - k_TL*TL = Ao*Ea
% Ia=Ea/Ra        Ia*kt - N*TL = 0       Vr = k_theta*theta

%k_TL=0 % No realimento en forma feedforward la estima de la perturbacion
k_TL=-(N/kt)*(k_Ia+Ao*Ra)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% TL como VE ==> TL[k+1]=TL[k]
% Modifico el modelo

A_mod=[1 0 0 0; Bd(:,2) Ad]
B_mod=[0;Bd(:,1)]
C_mod=[0 1 0 0;0 0 1 0;0 0 0 1]
D_mod=[0;0;0]

% Transformacion lineal T
T=eye(4)
At=T*A_mod*T^(-1);
Bt=T*B_mod;
Ct=C_mod*T^(-1);
Dt=D_mod;

At11=At(1,1)
At12=At(1,2:4)
At21=At(2:4,1)
At22=At(2:4,2:4)

Bt1=Bt(1)
Bt2=Bt(2:4)

% Calculo del estimador

polo_est=0.3;

% Aca surgen varias posibilidades porque hay una ecuacion y 3 incognitas:
%
%     1 - ( h1*At21(1) + h2*At21(2) + h3*At21(3) ) = polo_est
%
h1=0;
h3=0;
h2=(1-polo_est)/At21(2);

H=[h1 h2 h3]

Ae=At11-H*At21
Be=Bt1-H*Bt2
He=At12-H*At22+Ae*H

% Ganancias para implementar el estimador en Simulink
He_Ia=He(1);
He_omega=He(2);
He_theta=He(3);

H_Ia=H(1);
H_omega=H(2);
H_theta=H(3);

% Ruido en las mediciones
noise_power_Ia=0*1e-6;
noise_power_omega=0*1e-6;
noise_power_theta=0*1e-7;

