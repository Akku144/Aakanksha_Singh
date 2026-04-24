function [dvTot, dtTot] = transfer_cost(oe1, oe2, tNow, model)
% Fast semi-analytical EOR-style transfer estimate (Huang et al. inspired).
% oe = [h_km, e, i_deg, RAAN_deg, argp_deg, M_deg]
% This implementation follows the paper structure:
% 1) Eq.(4)-style transfer time from RAAN drift
% 2) Eq.(16)-style linear solution for [Delta a, Delta i, Delta Omega]
% 3) Eq.(17)-style two-impulse Delta-v estimate

mu = model.mu; Re = model.Re; J2 = model.J2;

% Numerical safeguards (defaults can be overridden by model fields).
if isfield(model, 'minDriftTimeSec')
    minDriftTimeSec = model.minDriftTimeSec;
else
    minDriftTimeSec = 86400; % 1 day floor to avoid phase-term blow-up
end
if isfield(model, 'maxDriftTimeSec')
    maxDriftTimeSec = model.maxDriftTimeSec;
else
    maxDriftTimeSec = inf; % optional practical cap
end
if isfield(model, 'kRangeEq4')
    kRangeEq4 = model.kRangeEq4;
else
    kRangeEq4 = -20:20; % Eq.(4) branch search over 2*pi*k
end
if isfield(model, 'maxPhaseFrac')
    maxPhaseFrac = model.maxPhaseFrac;
else
    maxPhaseFrac = 0.20; % cap |dAphase/a| contribution
end
if isfield(model, 'maxDvKmps')
    maxDvKmps = model.maxDvKmps;
else
    maxDvKmps = 20; % guard against nonphysical outliers
end

a1 = Re + oe1(1); e1 = oe1(2); i1 = deg2rad(oe1(3));
O1 = deg2rad(oe1(4)); w1 = deg2rad(oe1(5)); M1 = deg2rad(oe1(6));
a2 = Re + oe2(1); e2 = oe2(2); i2 = deg2rad(oe2(3));
O2 = deg2rad(oe2(4)); w2 = deg2rad(oe2(5)); M2 = deg2rad(oe2(6));

n1 = sqrt(mu / a1^3);
n2 = sqrt(mu / a2^3);
p1 = a1 * (1 - e1^2);
p2 = a2 * (1 - e2^2);

Odot1 = -1.5 * J2 * n1 * (Re / p1)^2 * cos(i1);
Odot2 = -1.5 * J2 * n2 * (Re / p2)^2 * cos(i2);
wDot1 = 0.75 * J2 * n1 * (Re / p1)^2 * (5 * cos(i1)^2 - 1);
wDot2 = 0.75 * J2 * n2 * (Re / p2)^2 * (5 * cos(i2)^2 - 1);
MDot1 = n1 + 0.75 * J2 * n1 * (Re / p1)^2 * sqrt(1 - e1^2) * (3 * cos(i1)^2 - 1);
MDot2 = n2 + 0.75 * J2 * n2 * (Re / p2)^2 * sqrt(1 - e2^2) * (3 * cos(i2)^2 - 1);

% Eq.(4)-style transfer duration driven by RAAN difference elimination.
% Search 2*pi*k branches and pick the smallest positive rendezvous time.
dO0 = wrapToPiLocal(O2 - O1);
dOdot = Odot2 - Odot1;
if abs(dOdot) < 1e-14
    dt = inf;
else
    dOmegaCandidates = dO0 + 2*pi*kRangeEq4(:);
    tCand = abs(dOmegaCandidates ./ dOdot);
    tCand = tCand(isfinite(tCand) & tCand > 0);
    if isempty(tCand)
        dt = inf;
    else
        dt = min(tCand);
    end
end
dt = min(dt, maxDriftTimeSec);
if ~isfinite(dt)
    dvTot = inf;
    dtTot = inf;
    return;
end
dtEff = max(dt, minDriftTimeSec);

% Differences at rendezvous epoch tf = tNow + dt without impulses.
O1f = O1 + Odot1 * (tNow + dt);
O2f = O2 + Odot2 * (tNow + dt);
dOmegaTf = wrapToPiLocal(O2f - O1f);

dA0 = a2 - a1;
dI0 = i2 - i1;
OmegaDot0 = Odot1;

% Eq.(16) linear solution for x1=Delta a / a, x2=Delta i, x3=Delta Omega.
si = sin(i1);
if abs(si) < 1e-6, si = sign(si + 1e-12) * 1e-6; end
den = (49 / 2) + (tan(i1)^2 / 2) + (2 / ((si^2) * (dtEff^2) * (OmegaDot0^2 + 1e-30)));
rhs = (dOmegaTf / dtEff) / (OmegaDot0 + 1e-30);
lam = rhs / den;
x1 = 7 * lam;
x2 = 0.5 * tan(i1) * lam;
x3 = -(2 / (si^2 * dtEff * (OmegaDot0 + 1e-30))) * lam;

dA = x1 * a1;
dI = x2;
dOmega = x3;
V = sqrt(mu / a1);

% Eq.(17)-style two-impulse estimate.
term0 = (dA / (2 * a1))^2 + dI^2 + ((dOmega * si) / 2)^2;
termf = ((dA - dA0) / (2 * a1))^2 + (dI - dI0)^2 + ((dOmega * si) / 2)^2;
dv0 = V * sqrt(max(term0, 0));
dvf = V * sqrt(max(termf, 0));

% Phase correction (Sec. III-B, Eq. 22/23 approximation).
u1f = w1 + M1 + (wDot1 + MDot1) * (tNow + dt);
u2f = w2 + M2 + (wDot2 + MDot2) * (tNow + dt);
dU = wrapToPiLocal(u2f - u1f);
dAphase = -2 * a1 * dU / max(3 * n1 * dtEff, 1e-12);
frac = abs(dAphase / a1);
if frac > maxPhaseFrac
    dAphase = sign(dAphase) * maxPhaseFrac * a1;
end
dvPhase = 0.5 * V * abs(dAphase / a1);

% Eccentricity-vector correction (Sec. III-C lightweight approximation).
w1f = w1 + wDot1 * (tNow + dt);
w2f = w2 + wDot2 * (tNow + dt);
dEx = e2 * cos(w2f) - e1 * cos(w1f);
dEy = e2 * sin(w2f) - e1 * sin(w1f);
dE = hypot(dEx, dEy);
dvEcc = V * dE;

dvTot = model.transferEta * (dv0 + dvf + dvPhase + dvEcc);
dtTot = dt;

if ~isfinite(dvTot) || dvTot < 0
    dvTot = inf;
end
if dvTot > maxDvKmps
    dvTot = inf;
end
end

function y = wrapToPiLocal(x)
y = mod(x + pi, 2 * pi) - pi;
end
