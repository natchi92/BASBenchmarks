% Component: Boiler
% author: Nathalie Cauchi
% -------------------------------------------------------
% \dot(Tsw,b) = (taub)^(-1)[-(Tsw,b(t)) + kb ]dt + sigmadW
% -------------------------------------------------------
% -------------------------------------------------------

function M = BoilerModel(Boiler,Ts)
% Defining corresponding matrices

Ac  = -1/Boiler.taub;
Bc  = 0;
Cc   = 1;
Fc  = Boiler.kb/Boiler.taub;

%%%%%%%%% Create stochastic model %%%%%%%%%%%%%%
M = createModel;
M = InitialiseModel(M,'l','s',Ac,Bc,Cc,Fc,[],[],Ts,Boiler.sigma);
M = createSymbModel(M);