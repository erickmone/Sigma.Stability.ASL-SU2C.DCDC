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

% kp=[0.036,  0.09,  0.15, 0.16,  0.19,  0.12,  0.076, 0.056, 0.057];
% ki=[270,    1.4e3, 460,  3.2e3, 5.2e3, 6.6e3, 6.5e3, 6.2e3, 5.4e3];

% kp=[0.036, 0.19, 0.057];
% ki=[270, 5.2e3, 5.4e3];

% kp=0.036;
% ki=270;

% kp=0.19;
% ki=5.2e3;
% 
kp=0.057;
ki=5.4e3;

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

figure(1)

tOut=0:0.00002:3e-3;
[splot,tOut]=step(Mc,tOut);
plot(tOut,splot)

% Tabla1Splot=[tOut,splot];
% save('Tabla3StepPlot.txt','Tabla1Splot', '-ASCII','-append');

% hold on
% K=vpa([kp ki],2)
% 
opts = bodeoptions('cstprefs');
opts.PhaseVisible = 'off';
opts.FreqUnits = 'Hz';
opts.XLabel.Interpreter = 'latex';
opts.YLabel.Interpreter = 'latex';
opts.Title.Interpreter = 'latex';
opts.Title.FontSize = 12;
opts.XLabel.FontSize = 12;
opts.YLabel.FontSize = 12;
opts.MagUnits = 'abs';
% 
figure(2)
% bodeplot(OL,opts)
% hold on

bodeplot(Mc,opts,'m')
win=linspace(0,2*pi*1e5,500);
%% por defaul mag=abs, wout=rad/sec 
[mag,phase,wout]=bode(Mc,win);
Tabla1Bode=[wout./(2*pi),mag(:)];
figure(3)
plot(wout./(2*pi),mag(:))
save('Tabla3Bode.txt','Tabla1Bode', '-ASCII','-append');
% 
