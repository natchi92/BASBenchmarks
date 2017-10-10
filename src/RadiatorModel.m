% Component: Radiator
% author: Nathalie Cauchi
% -------------------------------------------------------
% \dot(Trw,r) = [k0*wrad(Trad(t)-Trw,r(t)) +
% k1(Tz(t) - Trw,r(t))]dt + sigmadW
% -------------------------------------------------------
% -------------------------------------------------------

function M = RadiatorModel(Radiator,Ts)
% Defining corresponding matrices

Ac = -Radiator.k1;
Bc = [Radiator.k0 0 Radiator.k1];
Nc = [0 -Radiator.k0 0];
Cc = 1;

%%%%%%%%% Create stochastic model %%%%%%%%%%%%%%
M = createModel;
M = InitialiseModel(M,'b','s',Ac,Bc,Cc,[],Nc,[],Ts,Radiator.rw.sigma);
M = createSymbModel(M);