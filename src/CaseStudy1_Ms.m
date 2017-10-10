% Case study 1: Stochastic Model 
% author: Nathalie Cauchi
% -------------------------------------------------------
% x_d[k+1] = Ax_d[k] + Bu[k] + Q_d + SigmaW[k]
% y_d[k]   = [1 0 0 0; 0 1 0 0]
% x_d = [T_z1 T_z1 T_rw,rad1 T_rw,rad2]^T
% u   = T_sa
% -------------------------------------------------------
% -------------------------------------------------------

clc; clear; 

% Load parameters needed to build model
BASParameters;
Ts = 15;                 % Sample time (minutes)
T  = 4*24*3;             % Simulation 

% Steady state values
Tswb  = Boiler.Tswbss;              
Tsp   = Zone1.Tsp;             
Twss  = Zone1.Twss;
Trwrss= Radiator.Trwrss;
Pout1 = Radiator.Zone1.Prad;    
Pout2 = Radiator.Zone2.Prad;    
m     = Zone1.m;
w     = Radiator.w_r;

Ac = zeros(4,4);
Ac(1,1) = -(1/(Zone1.Rn*Zone1.Cz))-((Pout1*Radiator.alpha2)/(Zone1.Cz)) - ((m*Materials.air.Cpa)/(Zone1.Cz));
Ac(1,3) = (Pout1*Radiator.alpha2)/(Zone1.Cz);
Ac(2,2) = -(1/(Zone2.Rn*Zone2.Cz))-(Pout2*Radiator.alpha2)/(Zone2.Cz) - (m*Materials.air.Cpa)/(Zone2.Cz);
Ac(2,4) = (Pout2*Radiator.alpha2)/(Zone2.Cz);
Ac(3,1) = (Radiator.k1);
Ac(3,3) = -(Radiator.k0*w)-Radiator.k1;
Ac(4,2) = (Radiator.k1);
Ac(4,4) = -(Radiator.k0*w) -Radiator.k1;

Bc = [(m*Materials.air.Cpa)/(Zone1.Cz) (m*Materials.air.Cpa)/(Zone2.Cz) 0 0]';

d =(Twss/(Zone1.Rn*Zone1.Cz))+ (Radiator.alpha1)/(Zone1.Cz);
g =((Twss-1)/(Zone2.Rn*Zone2.Cz))+ (Radiator.alpha1)/(Zone1.Cz);
k = (Radiator.k0*w*Tswb);
p = (Radiator.k0*w*Tswb);

Fc = [d g k p]';

Cc = diag([1 1 0 0]);

Sigma =diag([Zone1.Tz.sigma Zone2.Tz.sigma Radiator.rw.sigma Radiator.rw.sigma]);
Z1m=createModel;

% Creation of symbolic deterministic model
Z1m=InitialiseModel(Z1m,'l','s',Ac,Bc,Cc,Fc,[],[],Ts,Sigma);
Z1m=createSymbModel(Z1m);

% Defining input signal 
tsa = 20;
Tsa = 18.*ones(T,1);
Tsa(32:48)=tsa.*ones(17,1);
Tsa(52:72)=tsa.*ones(21,1);
Tsa(32+96:48+96)=tsa.*ones(17,1);
Tsa(52+96:72+96)=tsa.*ones(21,1);
Tsa(32+96*2:48+96*2)=tsa.*ones(17,1);
Tsa(52+96*2:72+96*2)=tsa.*ones(21,1);
Tsa(32+96*3:48+96*3)=tsa.*ones(17,1);
Tsa(52+96*3:72+96*3)=tsa.*ones(21,1);

D   = ones(T,1);    
   
dWz        = sqrt(Z1m.dt).*randn(T,size(Z1m.A,2));         % Brownian increments
Z1m.dW     = dWz;

% Simulate model over given time horizon T 
Tz1_y      = runModel(Z1m,[20 20 Trwrss Trwrss]',Tsa,D,T);

% Plot Results
title={{'Zone 1 temperature (^oC)'},{'Zone 2 temperature (^oC)'}};
plotFigures(1:T,Tz1_y(1:2,1:T),title); 
