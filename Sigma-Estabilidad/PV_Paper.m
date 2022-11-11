%% Modeling of the Photovoltaic array

% The current-voltage characteristic of the PV array can be written as 

Ns = 1; % Number of series cells
Np = 1; % Number of parallel cells

T = 30 + 273.15; % PV array temperature in Kelvin 
q = 1.60217662e-19; % Charge of an electron C
K0 = 1.3805e-23; % Boltzmanns constant J/K
A = 1.12; % Ideality Factor of thge P-N junction [1,5]

Tr = 298; % Reference Temperature K
ior = 5.98e-8; % Reverse saturation current A
Lambda = 900; % Irradiance W/m^2
Ego = 1.2; % Band-Gap energy of the semiconductor eV
iscr = 3.45; % Short-circuit current reference A
Kl = 12e-4; % Temperature coefficient A/K

irs=ior*((T/Tr)^3)*exp(((1/Tr)-(1/T))*((q*Ego)/(K0*A)));
iph=(iscr+Kl*(T-Tr))*Lambda/1000;
ip = Np*(iph-irs*(exp((q*Vp)/(Ns*A*K0*T)-1)));

