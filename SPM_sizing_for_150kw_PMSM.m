%%%%%%%%%%%%%%%% Analytical Design of Surface-Mount PMSM 150KW %%%%%%%%%%%%%%%%%%
clear all; clc
%%%%%%%%%%%%% H64AEM %%%%%%%%%% Advanced Electrical Machines %%%%%%%%%%%%%%
%% Preliminary design of an SPM motor with inner rotor
% Implementation of the equations derived during the lectures

%%%%%%% Data given %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Torque =  150.78;             % Required rated torque [Nm]
Speed =  9500 * pi/30;     % Required rated mechanical pulsation [rad/sec]
n = 9500   %  Base Speed [rpm]
max_speed = 12000  %maximum speed in [RPM]
V = 400/sqrt(3) ;          % Maximum RMS phase voltage

%%%%%%% Data chosen or Assumed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = 72;                    % Number of slots 
m = 3 ;                    % Number of phases                      
p = 4;                     % Pole pairs 
fs = p* Speed /(2*pi);     % Electrical frequency at rated speed 
q = Z/(2*m*p)              % Number of slots per pole/phase    
Airgap = 1.2e-3;          % Airgap length [m]
Ka = 0.83;                 % Magnet span to pole span ratio
Bmg1 = 0.8;                 % B*sqrt(2)(Fundamental peak airgap flux density [T])
B_g  = Bmg1 /(4/pi * sin(Ka*pi/2))     % Peak airgap flux density [T]
B = Bmg1/sqrt(2)         % Fundamental RMS airgap flux density [T]
B_mag = 1.04;             % Magnets remanent flux density [T]
Kfill = 0.45;             % Slot Fill factor (area of copper/Area of slot)
Rho = 2.7e-8;            % Copper resistivity at 100C
B_yoke = 1.7;           % Saturation flux densities in the stator yoke (back iron) 
B_tooth = 1.68;          % Saturation flux densities in the stator tooth
PF = 0.88 ;             % Estimate of power factor 
J = 11.5e6 ;            % Current density in the slot {4e6 [A/m^2] = 4 [A/mm^2]}
H = 42.0e3 ;            % Electrical loading [kA/m] also defined as A in the lectures
AR = 1;                 % Rotor Aspect Ratio AR = L/D   
mu_m=1.02;              % Magnet relative permeability (1)

%%%%%%%%%% Magnet thickness calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lm = Airgap * mu_m*((B_g./B_mag) /(1 - (B_g./B_mag)))
Lm1 = Airgap * mu_m / ((B_mag/B_g)-1)           % Simpler way to express the PM lenght Lm [m]:

%%%%%%%%%% Fundamental Winding distribution factor calculation %%%%%%%%%%%%
Kw1 = sin(pi/6)/(q*sin(pi/(6*q)))

%%%%%%%%%% Stator, Rotor and slots dimensions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%D = ((2*Torque) / (H*B*AR*Kw1*pi))^(1/3); Lz = D*AR ; 
Lz = 0.15;                      % Active axial lenght [m]
omega = 2*pi*n/60               % Mechanical pulsation [rad/s]
P = Torque * omega              % Mechanical power expressed at the shaft [W]

% D_ = ((2*Torque) / (H*B*Kw1*pi*Lz))^(1/2)
% Rotor Diameter determined from the Torque equation expressed in [m]
D = sqrt(((2*Torque)/(H*B*Kw1*pi*Lz)))

% Diameter (rotor diameter) calculated based on the demonstration
% presented during the lectures

D_validation = sqrt((120*P)/(H*Bmg1*Lz*n*sqrt(2)*(pi^2)*Kw1))

% % Airgap empirical equation g is the airgap expressed in [mm]
% g = 6*(D)/(sqrt(2*p)*1000)

% Preliminary definition of the stator geometry dimensions:
peripheral_Speed = (D/2)*2*pi*max_speed/60  %peripheral speed in meters per sec [ms-1]

fi_p = D*Lz*Bmg1/p                     %Flux per pole in Wb

tau_s = pi * D / Z                  % Slot pitch
Wt = tau_s * Bmg1 / B_tooth         % Tooth width
tau_p = pi * D /(2*p)               % Pole pitch
Wbi = tau_p * Bmg1 /(2*B_yoke)      % Back iron thickness (more conservative)
                                   %(NOTE the Back Iron is 
							       % depending on the pole number
								   %(in other words depending
                                   % on the pole pitch))
                                        									
% Wc = fi_p/(2*B_yoke*Lz)            %Alternative equation for back iron thickness                                      
                                        
Scu_all = pi*D*H / (J*Kfill);     % Copper section all slots
S_cu = Scu_all/Z;                 % Copper section one slot
Dbi = sqrt(8*Scu_all/pi + D^2)     % Diameter at the back iron (slot bottom)
Dext = Dbi + 2*Wbi                 % External diameter (Normally called
                                    % Outer Diameter OD or SOD Stator Outer Diameter)

% Calculation of the number of turns

Lew=2.5*(D/p);                      % End winding lenght (classical empirical equation)
L_tot= Lz + Lew;                    % Active length including an estimated value
                                    % of end-winding length 

% LOSSES ESTIMATION
% Copper Losses
P_copper = Rho*L_tot* Scu_all*Kfill*J^2 ;    % Copper loss 
Power = Torque*omega;                      % Mechanical output power 

% With the strong hypothesis that Pjoule = Piron
I = (Power + 2*P_copper) / (3*V*PF);       % RMS phase current  
Eff = Power / (Power + 2*P_copper);        % Efficiency
Nph = H * pi *D /(2*m*I);                  % Number of turns per phase
N = Nph/(p*q);                             % Number of turns per slot
Nph=int32(Nph);  N=int32(N);               % Get integer numbers of turns 

%%%%%%%% RMS Back EMF %%%%%%%%%%%%%%%%%
E1 = 2*pi*fs/p * Kw1*Nph*B*D*Lz; 

%%Temerature rise%% 
A_ext = pi*Dext*L_tot; %+ pi*Dext^2/4 ;     % External area
h = (2*P_copper)/(70 * A_ext) ;              % Heat convection transfer coeff 
DT = (2*P_copper)/(h*A_ext);                  % Temperature differential 

%%%%%%%Spread sheet display %%%%%%%%%%%%%%%%%%%%%%%%
disp(['Inner stator Diameter: ',num2str(D*1e3),' mm']);
disp(['Active axial length: ',num2str(Lz*1e3),' mm']);
disp(['Total axial length (with EW): ',num2str(L_tot*1e3),' mm']);
disp(['Magnet thickness: ',num2str(Lm*1e3),' mm']);
disp(['External Diameter: ',num2str(Dext*1e3),' mm']);
disp(['Tooth width: ',num2str(Wt * 1000),' mm']);
disp(['Slot Depth: ',num2str((Dbi-D)/2 * 1000),' mm']);
disp(['Phase current: ',num2str(I),' A']);
disp(['Number of turns per phase: ',num2str(Nph),' turns']);
disp(['Number of turns per slot: ',num2str(N),' turns']);
disp(['EMF ',num2str(E1),' V']);
disp(['Efficiency ',num2str(Eff*100),' %']);
disp(['Temperature rise ',num2str(DT),' deg C']);