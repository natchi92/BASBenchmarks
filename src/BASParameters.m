% BAS Parameters used through-out benchmarks
% Each BAS component is defined as a structure containing
% all the associated parameters
% author: Nathalie Cauchi
% -------------------------------------------------------


%--------- Material properties   -------%

Materials.air.Cpa   = 1204/1e3;                                             % thermal capacity of air [J/kgK]
Materials.air.rhoa  = 1.2;                                                  % density of air [kg/m3]
Materials.water.Cpw = 4180/1e3;                                             % thermal capacity of water [J/kgK]
Materials.water.rhow= 1000/1e3;                                             % density of water [kg/m3]
Materials.concrete.Cn = 880/1e3;                                            % thermal capacity of concrete [J/kgK]
Materials.concrete.rhon = 2.4e3/1e3;                                        % density of concrete [kg/m3]


%--------- Boiler Parameters  -------%

% Gas Boiler, AMBI simulator

Boiler.taub  = 60*5;                                                         % Time constant of boiler [s]
Boiler.kb    = 348.15 - 273.15;                                              % Steady-state temperature of the boiler (75degC) [K]
Boiler.Tswbss= 75;
Boiler.sigma = 0.5;



Splitter.uv = 0;                                                            % Water flow Splitter [-]



%--------- Fan Coil Unit (FCU) Parameters -------%
FCU.Vfcu   = 0.06;                                                          % Volume of fcu [m3]
FCU.Afcu   = 0.26;                                                          % Contact area for conduction in fcu [m2]
FCU.mfcu_l = 0.16;                                                          % FCU mass air flow when fan is in low mode [m3/s]
FCU.mfcu_m = 0.175;                                                         % FCU mass air flow when fan is in med mode [m3/s]
FCU.mfcu_h = 0.19;                                                          % FCU mass air flow when fan is in high mode [m3/s]

%--------- Air Hanlding Unit (AHU) Parameters     -------%

AHU.sa.k0    = 2.0167e-05;                                                  % Constant k0 = Cpa/(Cpa*rho_a*Vahu)
AHU.sa.k1    = 0.0183;                                                      % Constant k1 = (UA)a/(Cpa*rho_a*Vahu)
AHU.rw.k0    = 0.0109;                                                      % Constant k0 = Cpw/(Cpw*rho_h*Vahu)
AHU.rw.k1    = 0.0011;                                                      % Constant k1 = (UA)a/(Cpw*rho_h*Vahu)
AHU.rw.alpha3= 0.0011;                                                      % Constant alpha_3
AHU.rw.Trwss = 35;                                                          % AHU return water steady state temperature [deg C]
AHU.w_a      = 1/60;                                                        % Nominal rate of water flow [m3/min]
AHU.w_max    = 10/60;                                                       % Max rate of water flow [m3/min]
AHU.sa.sigma =0.1;
AHU.rw.sigma =0.1;

Mixer.um     = 0.5;
%--------- Radiator Parameters-------%

Radiator.k0         = 0.0193;                                                 % Constant k0 = Cpw/(Cpw*rho_h*Vr)
Radiator.k1         = 8.900e-04;                                              % Constant k1 = (UA)r/(Cpw*rho_h*Vr)
Radiator.w_r        = 5/60;                                                   % Nominal rate of water flow [m3/min]
Radiator.w_max      = 10/60;                                                  % Max rate of water flow [m3/min]
Radiator.Zone1.Prad = 800/1000;
Radiator.Zone2.Prad = 600/1000;                                               % Rated output power of radiator [kW]
Radiator.alpha2     = 0.0250;                                                 % Coefficients of Qrad [-]
Radiator.alpha1     = -0.02399;
Radiator.Trwrss     = 35;                                                     % Radiator steady state temperature [deg C]
Radiator.rw.sigma   = 0.1;

%--------- Zone Parameters    -------%
%--------- Zone 1             -------%

Zone1.Cz   =  51.5203;                                                       % Thermal capacitance of zone [J/kgK]
Zone1.Cn   =  52.3759;                                                       % Thermal capacitance of neighbouring wall [J/kgK]
Zone1.Rout =  7.20546;                                                       % Resistance of walls connected to outside [K/W]
Zone1.Rn   =  5.4176;                                                        % Resistance of walls connected to neighbouring wall [K/W]
Zone1.mu   =  0.000199702104146;                                             % Q_occ coefficient in \muCO_2 + \zeta [-]
Zone1.zeta =  0.188624079752965;
Zone1.alpha=  0.044;                                                         % Absorptivity coefficient of wall [W/m2K]
Zone1.A_w  =  1.352;                                                         % Area of window [m2]
Zone1.iota =  0.359444194048431;                                             % Q_solar coefficient in \alpha A_w(\iota T_out + \gamma) [-]
Zone1.gamma=  1.622446039686184;
Zone1.m    = 10/60;                                                          % Fixed mass supply air [m3/min]
Zone1.Twss = 18;                                                             % Zone 1 walls steady-state temperature [deg C]
Zone1.Tsp  = 20;                                                             % Zone 1 desired temperature set-point  [deg C]
Zone1.Tz.sigma = 0.02;
Zone1.Te.sigma = 0.01;

%--------- Zone 2             -------%
Zone2.Cz   =  50.2437;                                                       % Thermal capacitance of zone [J/kgK]
Zone2.Cn   =  53.3759;                                                       % Thermal capacitance of neighbouring wall [J/kgK]
Zone2.Rout =  7.4513;                                                        % Resistance of walls connected to outside [K/W]
Zone2.Rn   =  5.7952;                                                        % Resistance of walls connected to neighbouring wall [K/W]
Zone2.mu   =  5.805064e-6;                                                   % Q_occ coefficient in \muCO_2 + \zeta [-]
Zone2.zeta =  -0.003990;
Zone2.alpha=  0.044;                                                         % Absorptivity coefficient of wall [W/m2K]
Zone2.A_w  =  10.797;                                                        % Area of window [m2]
Zone2.iota =  0.03572;                                                       % Q_solar coefficient in \alpha A_w(\iota T_out + \gamma) [-]
Zone2.gamma=  0.06048;
Zone2.m    =  10/60;                                                         % Fixed mass supply air [m3/min]
Zone2.Twss =  18;                                                            % Zone 2 walls steady-state temperature [deg C]
Zone2.Tsp  =  20;                                                            % Zone 2 desired temperature set-point  [deg C]



Zone2.Tz.sigma = .02;

Zone2.Te.sigma = .01;
