% Case study 2: Model with large number of continuous variables
% author: Nathalie Cauchi
% -------------------------------------------------------
% x_c[k+1] = A_cx_c[k] + B_cu_c[k] +F_cd_c[k] + Q_c
% y_c[k]   = [1]
% x_c = [T_z1 ]^T
% u_c = T_sa
% d_c =[T_out CO2_1 T_rw,r1]^T
% -------------------------------------------------------
% -------------------------------------------------------

clc; clear; 

% Load parameters needed to build model
BASParameters;
Ts = 15;                 % Sample time (minutes)
T  = 4*24*3;             % Simulation 

% Steady state values
Tsp   = Zone1.Tsp;             
Trwass= AHU.rw.Trwss;
Pout1 = Radiator.Zone1.Prad;    
m1    = Zone1.m;
w     = Radiator.w_r;

% Defining Deterministic Model corresponding matrices
Ac      = zeros(1);
Ac(1,1) = -(1/(Zone1.Rout*Zone1.Cn))-Materials.air.Cpa*m1/Zone1.Cz - (Pout1*Radiator.alpha2)/Zone1.Cz;

Bc      = zeros(1,1);
Bc(1,1) = m1*Materials.air.Cpa/Zone1.Cz;

Fc = zeros(1,4);
 
Fc(1,2) = Zone1.mu/Zone1.Cz;
Fc(1,3) = Pout1*Radiator.alpha2/Zone1.Cz;
Fc(1,4) = (Pout1*Radiator.alpha1 + Zone1.zeta)/Zone1.Cz + Zone1.alpha*Zone1.A_w*Zone1.gamma/(Zone1.Cn) + (AHU.rw.alpha3*Trwass)/Zone1.Cn;

Cc = ([1 ]);
 
% Creation of symbolic deterministic model
Z1m=createModel;
Z1m=InitialiseModel(Z1m,'l','d',Ac,Bc,Cc,Fc,[],[],Ts,0);
Z1m=createSymbModel(Z1m);

% Defining input signal 
Tsa                 =18.*ones(T,1);
Tsa(32:48)          =20.*ones(17,1);
Tsa(52:72)          =20.*ones(21,1);
Tsa(32+96:48+96)    =20.*ones(17,1);
Tsa(52+96:72+96)    =20.*ones(21,1);
Tsa(32+96*2:48+96*2)=20.*ones(17,1);
Tsa(52+96*2:72+96*2)=20.*ones(21,1);
Tsa(32+96*3:48+96*3)=20.*ones(17,1);
Tsa(52+96*3:72+96*3)=20.*ones(21,1);

% Defining disturbance signal
Tout = 1*randn(T,1) + 9;
CO2_1= 100*randn(T,1) + 500;
Trwr1= 5*randn(T,1) + 35;

D = [Tout CO2_1 Trwr1 ones(T,1)];

% Simulate model over given time horizon T 
Tz1_y      = runModel(Z1m,[ 20 ]', Tsa,D,T);

% Plot Results
title={{'Zone 1 temperature (^oC)'},{'Input supply air (^oC)'}};
plotFigures(1:T,[Tz1_y(1,1:T); Tsa(1:T)' ],title); 
