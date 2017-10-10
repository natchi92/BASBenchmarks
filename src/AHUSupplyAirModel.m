% Component: AHU heating coil
% author: Nathalie Cauchi
% -------------------------------------------------------
% \dot(Tsa,a) = (k0*m(t)*(Td(t) - Tsa(t)) + k1(Tz(t) - Tsa(t)))dt
% + sigmadW
% -------------------------------------------------------
% -------------------------------------------------------

function M = AHUSupplyAirModel(AHU,Ts)
% Defining corresponding matrices
AHU.sa.k1=AHU.sa.k1/60;
Ac = -AHU.sa.k1;
Bc = [AHU.sa.k0 0 AHU.sa.k1];
Nc = [0 -AHU.sa.k0 0];
Cc = 1;


%%%%%%%%% Create stochastic model %%%%%%%%%%%%%%
M = createModel;
M = InitialiseModel(M,'b','s',Ac,Bc,Cc,[],Nc,[],Ts,AHU.sa.sigma);
M = createSymbModel(M);