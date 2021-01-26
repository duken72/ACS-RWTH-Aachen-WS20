% Clear all previous parameters and close all figures
clear all;
close all;
clc;

% General conversion factors
HG2PA        = 133.322368;       % mmHg => Pa
PA2HG        = 1.0/HG2PA;        % Pa => mmHg
LMIN2M3      = 1/(1.0E03*60);    % l/min => m^3  
M32LMIN      = 1.0E03*60;        %  m^3  => l/min  
SAMPLE       = 0.001;            % Sampling time [s]
RADS2TURNMIN = 60/(2*pi);        % rad/s => RPM

% Parameter for the gear pump subsystem
Lgp     = 1.44*(10^-3); % Motor inductance [H]
Rgp     = 0.48; % Motor electrical resistance [Ohm]
Kgpemf  = 72.4*(10^-3); % Electro-magnetic-mechanical coupling [Vs]
w0      = 293.2153; % Rotational neutral gear speed [rad/s]
I0      = 1.15; % Current neutral gear [A]
MR1     = Kgpemf*I0/w0; % Viscous damping [VAs^2]
Jgp     = 0.293*(10^-3); % Inertia gear pump [kgm^2]
lgp     = 0.025; % Gear length [m]
ri      = 0.012; % Inner gear radius [m]
ra      = 0.017; % Outer gear radius [m]
Mlgp    = lgp*(ra^2 - ri^2); % Load moment constant hydraulic coupl. [m^3]
Imaxzp  = 7.5; % Maximum current gear pump [A]
Rhyd    = 572.46e6;         % Hydraulic resistance tubing [Pa/(m^3/s)]
Lhyd    = 14.92e6;          % Inertance tubing [kg/m^4]
Vgpsv   = 1.59*(10^-6); % Swallowing capacity of gear pump [m^3]
tr0     = 44.9*(10^-6); %[N]
tr1     = 0.3*(10^-3); %[Ns/m]

% Parameter for the voice coil actuator subsystem
Lvca    = 2.9*1.0E-03; % inductance [H]
Rvca    = 2.6; % electrical resistance [Ohm]
Kvcaemf = 35.14; % Electromechanic-coupling [N/A]
Kemfvca = 35.14; % Electromagnetic feedback [V*s/m]
mbel    = 1.1; % bellows mass [kg]
cbel    = 741.22; % bellows spring constant [N/m]
dbel    = 1.8; % viscous damping [Ns/m]
Abel    = 0.0047; % new bellow surface [m^2]

% Pressure compartment parameters
dcomp   = 0.08; % Diameter of pressure chamber [m]
lcomp   = 0.06; % Length of pressure chamber [m]
lbel    = 0.1; % Length of bellow max. stretch [m]
Acomp   = dcomp^2*pi/4;
V0comp  = Acomp*lcomp; % Volume of pressure chamber [m^3]
V0bel   = Abel*lbel/2; % Volume of bellow [m^3]
V0      = V0comp + V0bel; % Total volume [m^3]
EH2O    = 2.15*1.0E09; % Bulk modulus of water [Pa]

% General parameters pre-tubing and flow measurement filter
R1      = 50.0/8.0*HG2PA/(0.001/60); % Hyd. resistance of left tubing [Pa/(m^3/s)]
RS      = 30.0/15.0*HG2PA/(0.001/60); % Hyd. resistance of left tubing [Pa/(m^3/s)]
Rbf     = 1000000000000; %[Pa/(m^3/s)]
rhoH2O  = 1.0E03; % Density of water [kg/m^3]
rhoB    = 1.055E03; % Desnity of blood [kg/m^3]
g       = 9.81; % Earth gravity [m/s^2]
l0      = 1.3; % Length of tubing to water reservoir [m]
ct      = 25e-12;       % Compliance of the tubing  [s^2*m^4/Kg]
cpc     = 2.15*(10^9); %bulk modulus of the MCL fluid including the compressibility of the pressure champer [Pa]
h       = 1.3;  %[m];
ps      = 14.4*(10^3);% [Pa]     

% Initial conditions
omega_zp_init = 0.0; % [rad/s]
im_zp_init    = 0.0; % [A]
phi_zp_init   = 0.0; % [rad]

v_vca_init    = 0.0; % [m/s]
im_vca_init   = 0.0; % [A]
x_vca_init    = 0.0; % [m] [0, 0.05]

p_init = 0.0; % Initial pressure pressure chamber [Pa]