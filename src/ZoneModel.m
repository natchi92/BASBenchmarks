% Component: Zone
% author: Nathalie Cauchi
% -------------------------------------------------------
% \dot(Tz) = (Cz)^(-1)[ (Te(t) - Tz(t))/ Rout + Pout(alpha2(Trad(t)-Tz(t)) +p10) + muCO2(t) + zeta + ma(t)Cpa(Tsa(t) - Tz(t))]
% \dot(Te) = (Cn1)^(-1) [(Tout(t) - Te(t)/Ro + alphaA1(a Tout(t) + b) + (Tz(t) - Te(t)/Riw]
% -------------------------------------------------------
% -------------------------------------------------------

function M = ZoneModel(Materials,Radiator,Zone,Ts)
% Defining corresponding matrices

Cz1 = Zone.Cz;
Ce1 = Zone.Cn;
Rout= Zone.Rout;
Rn  = Zone.Rn;

a   = -(1/(Rn*Cz1))- (Zone.m*Materials.air.Cpa/Cz1) - (Radiator.Zone1.Prad*Radiator.alpha2/Cz1);
b   = 1/(Rn*Cz1);
c   = 1/(Rn*Ce1);
d   = -1/(Rout*Ce1) -  1/(Rn*Ce1);

Ac  = [a b; c d];


%% u = Tsa, Trad, TSP
e  = (Zone.m*Materials.air.Cpa/Cz1);
f  = (Radiator.Zone1.Prad*Radiator.alpha2/Cz1);

Bc =[e f; 0 0];


%% d = CO2, Tout, K
k = Zone.mu/ Cz1;
m = (Radiator.Zone1.Prad*Radiator.alpha1 + Zone.zeta)/Cz1;
p = (Zone.alpha*Zone.A_w*Zone.iota/Ce1) + (1/(Rout*Ce1));
q = (Zone.alpha*Zone.A_w*Zone.gamma/Ce1);

Fc = [k 0 m ; 0 p q];

Cc = [1 0;0 1];

%%%%%%%%% Create Models %%%%%%%%%%%%%%
M    = createModel;

%%%%%%%%% Create stochastic model %%%%%%%%%%%%%%
M = InitialiseModel(M,'l','s',Ac,Bc(1:2,1:2),Cc,Fc,[],[],Ts,[Zone.Tz.sigma Zone.Te.sigma]);
M = createSymbModel(M);

