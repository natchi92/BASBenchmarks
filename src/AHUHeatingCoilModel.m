% Component: AHU heating coil
% author: Nathalie Cauchi
% -------------------------------------------------------
% \dot(Trw,a) = (k0*wahu(t)*(Tsw,b(t) - Trw,a(t)) + k1(Tz(t) - Trw,a(t)))dt
% + sigmadW
% -------------------------------------------------------
% -------------------------------------------------------

function M = AHUHeatingCoilModel(AHU,Ts)
% Defining corresponding matrices

Ac = -AHU.rw.k1;
Bc = [AHU.rw.k0 0 AHU.rw.k1];
Nc = [0 -AHU.rw.k0 0];
Cc = 1;


%%%%%%%%% Create stochastic model %%%%%%%%%%%%%%
M = createModel;
M = InitialiseModel(M,'b','s',Ac,Bc,Cc,[],Nc,[],Ts,AHU.rw.sigma);
M = createSymbModel(M);