% Benchmark: Stochastic Models representing the whole BAS set-up being
% taken under consideration (Fig.1)
% Components: 1 Boiler, 1 AHU, 1 mixer, 1 collector, 2 zones, 1 radiator in
% each zone
% author: Nathalie Cauchi
% -------------------------------------------------------------------------
clc; clear;
% Load parameters needed to build model
BASParameters;
Ts = 15;                 % Sample time (minutes)
T  = 4*24*3;             % Simulation
num= 24*4;               % number of minutes in 1 day
Tsp= 20;                 % Zone set-point
Tmax=21;                 % Max temperature for comfort
Tmin=19;                 % Min temperature for comfort
% Definining models of the individual components making up the BAS setup
Tswbm = BoilerModel(Boiler,Ts);                         % Boiler
Trwrm = RadiatorModel(Radiator,Ts);                     % Radiator 1
Trwrm2= RadiatorModel(Radiator,Ts);                     % Radiator 2
Trwm  = AHUHeatingCoilModel(AHU,Ts);                    % AHU heating coil
Tsam  = AHUSupplyAirModel(AHU,Ts);                      % AHU supply air
Tz1m  = ZoneModel(Materials,Radiator,Zone1,Ts);         % Zone 1
Tz2m  = ZoneModel(Materials,Radiator,Zone2,Ts);         % Zone 2
Wm    = ValveModel(Radiator.w_max);                     % Valve (radiator)
Mm    = MixerCollectorModel(0.5,2);                     % Mixer
Cm    = MixerCollectorModel(0.75,2);                       % Collector

% Initialise Model parameters
Tswb = zeros(1,T);
Trwr = zeros(1,T);
Trwr2= zeros(1,T);
Trwa = zeros(1,T);
Trwb = zeros(1,T);
Tsa  = zeros(1,T);
Tz1  = zeros(2,T);
Tz2  = zeros(2,T);
Td   = zeros(1,T);

% Initial values
Tswb(1,1)   =35;
Trwr(1,1)   =30;
Trwr2(1,1)  =30;
Trwa(1,1)   =30;
Tsa(1,1)    =19;
Tz1(1:2,1)  =[15 13]';
Tz2(1:2,1)  =[18 13]';

% Defining input signals
wsw_rad1             = (8/60)*ones(1,T);
wsa                  = (8/60)*ones(1,T);
wsw_AHU              = wsw_rad1;
wsw_rad2             = wsw_rad1;

CO2_1                = abs(repmat(500,1,T) + 100*randn(1,T));
CO2_2                = abs(repmat(500,1,T) + 100*randn(1,T));
Tout                 = repmat(9,1,T) + 1*randn(1,T);
Td(1,1)              = Mm(Tout(1,1),[Tz1(1,1) Tz2(1,1)]);
Trwr(1,1)            = Cm(Trwa(1,1), [Trwr(1,1) Trwr2(1,1)]);

% Simulation of whole BAS setup with basic condition based control on the
% water or air flow rates

% Models are stochastic
rng(100);
dWb   = sqrt(Tswbm.dt)*randn(size(Tswbm.A,2),T);         % Brownian increments
dWrwr = sqrt(Trwrm.dt)*randn(size(Trwrm.A,2),T);         % Brownian increments
dWrwr2= sqrt(Trwrm2.dt)*randn(size(Trwrm2.A,2),T);       % Brownian increments
dWrwa = sqrt(Trwm.dt)*randn(size(Trwm.A,2),T);           % Brownian increments
dWsa  = sqrt(Tsam.dt)*randn(size(Tsam.A,2),T);           % Brownian increments
dWz1  = sqrt(Tz1m.dt)*randn(size(Tz1m.A,2),T);           % Brownian increments
dWz2  = sqrt(Tz2m.dt)*randn(size(Tz2m.A,2),T);           % Brownian increments

for i=2:T
    Tswbm.dW  = dWb(1,i);
    Trwrm.dW  = dWrwr(1,i);
    Trwrm2.dW = dWrwr(1,i);
    Trwm.dW   = dWrwa(1,i);
    Tsam.dW   = dWsa(1,i);
    Tz1m.dW   = dWz1(1,i);
    Tz2m.dW   = dWz2(1,i);
    
    % Simple condition based control
    if Tz1(1,i-1) >= Tmax && Tz2(1,i-1) >=Tmax
        wsw_rad1(1,i) = Wm(0);
        wsw_rad2(1,i) = Wm(0);
        wsw_AHU(1,i)  = 0;
        wsa(1,i)      = 0;
        Tswb(1,i)     = 0;
    elseif Tz1(1,i-1) >= Tmax && Tz2(1,i-1)< Tmax && Tz2(1,i-1)>Tmin
        wsw_rad1(1,i) = Wm(0.5);
        wsw_rad2(1,i) = Wm(0.5);
        wsa(1,i)  = 0;
        Tswb(1,i)  = runModel(Tswbm,Tswb(1,i-1),0,1,1); % Boiler SW
    elseif Tz2(1,i-1) >= Tmax && Tz1(1,i-1)< Tmax && Tz1(1,i-1)>Tmin
        wsw_rad1(1,i) = Wm(0.5);
        wsw_rad2(1,i) = Wm(0.5);
        wsa(1,i)      = 0;
        Tswb(1,i)  = runModel(Tswbm,Tswb(1,i-1),0,1,1); % Boiler SW
        wsa(1,i)      = 10/60;
    elseif (Tz1(1,i-1) <= Tsp || Tz2(1,i-1) <= Tsp)
        wsw_AHU(1,i)  = 10/60;
        wsw_rad1(1,i) = Wm(1);
        wsw_rad2(1,i) = Wm(1);
        Tswb(1,i)  = runModel(Tswbm,Tswb(1,i-1),0,1,1); % Boiler SW
    end
    
     
    
    Td(1,i-1 ) = Mm(Tout(1,i-1), [Tz1(1,i-1) Tz2(1,i-1)]); % Collector
    
    Tsa(1,i)   = runModel(Tsam,Tsa(1,i-1),[(wsa(1,i-1).*(Td(1,i-1))) wsa(1,i-1) Tz1(1,i-1)],[],1);
    Trwr(1,i)  = runModel(Trwrm,Trwr(1,i-1),[wsw_rad1(1,i-1).*Tswb(1,i), wsw_rad1(1,i-1), Tz1(1,i-1)],[],1); % RAD 1
    Trwr2(1,i) = runModel(Trwrm2,Trwr2(1,i-1),[wsw_rad2(1,i-1).*Tswb(1,i), wsw_rad2(1,i-1), Tz2(1,i-1)],[],1); % RAD 2
    Trwa(1,i)  = runModel(Trwm,Trwa(1,i-1),[wsw_AHU(1,i-1).*Tswb(1,i), wsw_AHU(1,i-1), (Tz1(1,i-1)+Tz2(1,i-1))/2],[],1); % AHU return water
    
    
    Trwb(1,i-1)= Cm(Trwa(1,i-1), [Trwr(1,i-1) Trwr2(1,i-1)]); % Collector
    Tz1(:,i)   = runModel(Tz1m,[Tz1(1,i-1) Tz1(2,i-1)]',[Tsa(1,i-1) Tswb(1,i-1) ],[CO2_1(1,i-1) Tout(1,i-1) 1],1);
    Tz2(:,i)   = runModel(Tz2m,[Tz2(1,i-1) Tz2(2,i-1)]',[Tsa(1,i-1) Tswb(1,i-1) ],[CO2_2(1,i-1) Tout(1,i-1) 1],1);
    
end

% Plotting simulation output of whole BAS set-up
n =1;
title={{'Boiler supply water temperature (^oC)'},{'Boiler return water temperature (^oC)'}};
plotFigures(1:n:T, [Tswb(1,1:n:T); Trwb(1,1:n:T)],title);
title={{'Return water temperature of radiator 1 (^oC)'},{'Return water temperature of radiator 2 (^oC)'}};
plotFigures(1:n:T,[Trwr(1,1:n:T); Trwr2(1,1:n:T)],title);
title = {{'Return water temperature of AHU (^oC)'},{'Supply air temperature of AHU (^oC)'}};
plotFigures(1:n:T,[Trwa(1,1:n:T); Tsa(1,1:n:T)],title);
title = {{'Zone 1 temperature (^oC)'},{'Zone 2 temperature (^oC)'}};
plotFigures(1:n:T,[Tz1(1,1:n:T); Tz2(1,1:n:T)],title);
title = {{'AHU valve position'}, {'Radiators valve position'}};
plotFigures(1:n:T,[wsw_AHU(1,1:n:T); wsw_rad1(1,1:n:T)],title);
title = {{'Outside Temperature'}, {'CO_2 levels'}};
plotFigures(1:n:T,[Tout(1,1:n:T); CO2_1(1,1:n:T)],title);
