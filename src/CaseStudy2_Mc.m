% Case study 2: Model with large number of continuous variables
% author: Nathalie Cauchi
% -------------------------------------------------------
% x_c[k+1] = A_cx_c[k] + B_cu_c[k] +F_cd_c[k] + Q_c
% y_c[k]   = [1 0 0 0 0 0 0]
% x_c = [T_z1 T_z2 T_w5 T_w6 T_w2 T_w3 T_w7]^T
% u_c = T_sa
% d_c =[T_out T_hall CO2_1 CO2_2 T_rw,r1 T_rw,r2]^T
% -------------------------------------------------------
% -------------------------------------------------------

%clc; clear; 

% Load parameters needed to build model
BASParameters;
Ts = 15;                 % Sample time (minutes)
T  = 4*24*3;             % Simulation 

% Steady state values
Tsp   = Zone1.Tsp;             
Trwass= AHU.rw.Trwss;
Pout1 = Radiator.Zone1.Prad;    
Pout2 = Radiator.Zone2.Prad;    
m1    = Zone1.m;
m2    = Zone2.m;
w     = Radiator.w_r;

% Defining Deterministic Model corresponding matrices
Ac      = zeros(7);
Ac(1,1) = -(3/(Zone1.Rn*Zone1.Cz))-Materials.air.Cpa*m1/Zone1.Cz - (Pout1*Radiator.alpha2)/Zone1.Cz;
Ac(1,3) = 1/(Zone1.Rn*Zone1.Cz);
Ac(1,5) = 1/(Zone1.Rn*Zone1.Cz);
Ac(1,7) = 1/(Zone1.Rn*Zone1.Cz);
Ac(2,2) = -(3/(Zone2.Rn*Zone2.Cz))-Materials.air.Cpa*m2/Zone2.Cz - (Pout2*Radiator.alpha2)/Zone2.Cz;
Ac(2,4) = 1/(Zone2.Rn*Zone2.Cz);
Ac(2,6) = 1/(Zone2.Rn*Zone2.Cz);
Ac(2,7) = 1/(Zone2.Rn*Zone2.Cz);
Ac(3,1) = 1/(Zone1.Rn*Zone1.Cn);
Ac(3,3) = -(1/(Zone1.Rn*Zone1.Cn))-(1/(Zone1.Rout*Zone1.Cn)) - (AHU.rw.alpha3)/Zone1.Cn;             
Ac(4,2) = 1/(Zone2.Rn*Zone2.Cn);
Ac(4,4) = -(1/(Zone2.Rn*Zone2.Cn))-(1/(Zone2.Rout*Zone2.Cn)) - (AHU.rw.alpha3)/Zone2.Cn;
Ac(5,1) = 1/(Zone1.Rn*Zone1.Cn);
Ac(5,5) = -(1/(Zone1.Rn*Zone1.Cn))-(1/(Zone1.Rout*Zone1.Cn)) - (AHU.rw.alpha3)/Zone1.Cn;
Ac(6,2) = 1/(Zone2.Rn*Zone2.Cn);
Ac(6,6) = -(1/(Zone2.Rn*Zone2.Cn))-(1/(Zone2.Rout*Zone2.Cn)) - (AHU.rw.alpha3)/Zone2.Cn;
Ac(7,1) = 1/(Zone1.Rn*Zone1.Cn);
Ac(7,2) = (1/(Zone2.Rn*Zone2.Cn));
Ac(7,7) = -(1/(Zone1.Rn*Zone1.Cn))-(1/(Zone2.Rn*Zone2.Cn)) - (AHU.rw.alpha3)/((Zone1.Cn+Zone2.Cn)/2);

Bc      = zeros(7,1);
Bc(1,1) = m1*Materials.air.Cpa/Zone1.Cz;
Bc(2,1) = m2*Materials.air.Cpa/Zone2.Cz;

Fc = zeros(7,7);
 
Fc(1,3) = Zone1.mu/Zone1.Cz;
Fc(1,5) = Pout1*Radiator.alpha2/Zone1.Cz;
Fc(1,7) = (Pout1*Radiator.alpha1 + Zone1.zeta)/Zone1.Cz;
Fc(2,4) = Zone2.mu/Zone2.Cz;
Fc(2,6) = Pout2*Radiator.alpha2/Zone2.Cz;
Fc(2,7) = (Pout2*Radiator.alpha1 + Zone2.zeta)/Zone2.Cz;
Fc(3,1) = Zone1.alpha*Zone1.A_w*Zone1.iota/(Zone1.Cn) + 1/(Zone1.Rout*Zone1.Cn);
Fc(3,7) = Zone1.alpha*Zone1.A_w*Zone1.gamma/(Zone1.Cn) + (AHU.rw.alpha3*Trwass)/Zone1.Cn;
Fc(4,1) = Zone2.alpha*Zone2.A_w*Zone2.iota/(Zone2.Cn) + 1/(Zone2.Rout*Zone2.Cn);
Fc(4,7) = Zone2.alpha*Zone2.A_w*Zone2.gamma/(Zone2.Cn) + (AHU.rw.alpha3*Trwass)/Zone2.Cn;
Fc(5,2) = 1/(Zone1.Rout*Zone1.Cn);
Fc(6,2) = 1/(Zone2.Rout*Zone2.Cn);
Fc(5,7) = (AHU.rw.alpha3*Trwass)/Zone1.Cn;
Fc(6,7) = (AHU.rw.alpha3*Trwass)/Zone2.Cn;
Fc(7,7) = (AHU.rw.alpha3*Trwass)/((Zone1.Cn +Zone2.Cn)/2);


Cc = ([1 0 0 0 0 0 0 ]);
 
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
Thall= 1*randn(T,1) + 15;
CO2_1= 100*randn(T,1) + 500;
CO2_2= 100*randn(T,1) + 500;
Trwr1= 5*randn(T,1) + 35;
Trwr2= 5*randn(T,1) + 35;

D = [Tout Thall CO2_1 CO2_2 Trwr1 Trwr2 ones(T,1)];

% Simulate model over given time horizon T 
Tz1_y      = runModel(Z1m,[ 20 20 18 18 18 18 18]', Tsa,D,T);

% Plot Results
title={{'Zone 1 temperature (^oC)'},{'Input supply air (^oC)'}};
plotFigures(1:T,[Tz1_y(1,1:T); Tsa(1:T)' ],title); 


