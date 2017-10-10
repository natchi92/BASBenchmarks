% Component: Collector / Mixer
% author: Nathalie Cauchi
% -------------------------------------------------------
% Collector:
% Trw,b(t) = uvTrw,a(t)  + (1-uv)*sum(Twr_i)/i
% Mixer: 
% Td(t) = udTout + (1-ud)*sum(Tz_i)/i
% -------------------------------------------------------
% -------------------------------------------------------

function M = MixerCollectorModel(Ratio,n)


M = @(T1,T2) Ratio*T1 + (1-Ratio)*sum(T2)/n;