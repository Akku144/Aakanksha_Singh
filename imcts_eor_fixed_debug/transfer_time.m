function dt = transfer_time(oeC, oeT, P)
% Ye Eq. (4): dt = (Omega_j - Omega_i)/(Omegadot_j - Omegadot_i)
rateC = rates_j2(oeC, P);
rateT = rates_j2(oeT, P);
dOm = wrapToPiLocal(oeT(4) - oeC(4));
den = rateT(1) - rateC(1);
if abs(den) < 1e-12
    dt = 5*P.day;
else
    dt = dOm / den;
end
% enforce positive practical transfer by adding relative RAAN-drift cycles
period = 2*pi / max(abs(den),1e-12);
while dt <= 0, dt = dt + period; end
dt = max(dt, 0.5*P.day);
end

function y = wrapToPiLocal(x)
y = mod(x+pi,2*pi)-pi;
end
