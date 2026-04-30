function dvopt = proxy_cost(oeC0, oeT0, dt, P, Nmax)
% Huang et al. EOR-style estimator.
% Implements the published sequence:
% 1) solve linear KKT Eq. (14)-(16)
% 2) solve nonlinear a-i-RAAN problem Eq. (19)/(20)
% 3) phase correction Eq. (21)-(24), scan N
% 4) eccentricity correction via direct constrained Eq. (30)
%
% Requires Optimization Toolbox for fmincon/fsolve.

if nargin < 5, Nmax = 5; end
if dt <= 0 || ~isfinite(dt), dvopt = inf; return; end

oeC = propagate_oe_j2(oeC0, dt, P);
oeT = propagate_oe_j2(oeT0, dt, P);

a=oeC0(1); inc=oeC0(3); V=sqrt(P.mu/a);
rates0 = rates_j2(oeC0, P); Om0dot = rates0(1);

dOm0 = wrapToPiLocal(oeT(4) - oeC(4));
da0 = oeT(1) - oeC0(1);
di0 = oeT(3) - oeC0(3);

% Eq. 14-16 linear KKT seed
den = 49/2 + tan(inc)^2/2 + 2/(sin(inc)^2*(dt*Om0dot)^2);
lam = (dOm0/(dt*Om0dot)) / den;
xseed = [7*lam; 0.5*tan(inc)*lam; -2*lam/(sin(inc)^2*dt*Om0dot)];

opts = optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',80, ...
                    'OptimalityTolerance',1e-9,'ConstraintTolerance',1e-9);

% Eq. 19 minimization
fun19 = @(x) eor_cost_aiOmega(x, a, inc, da0, di0, V);
nonl19 = @(x) con_aiOmega(x, oeC0, dOm0, dt, P);
lb = [(P.amin-a)/a; -0.2; -2*pi]; ub = [0.2; 0.2; 2*pi];
try
    xstar = fmincon(fun19, xseed, [],[],[],[], lb, ub, nonl19, opts);
catch
    xstar = xseed;
end

dvbest = inf;
for N = -Nmax:Nmax
    xfix1 = phase_fixed_x1(xstar, oeC0, oeT, dt, P, N);
    if xfix1 < (P.amin-a)/a || xfix1 > 0.2, continue; end
    fun24 = @(y) eor_cost_aiOmega([xfix1;y(1);y(2)], a, inc, da0, di0, V);
    nonl24 = @(y) con_aiOmega([xfix1;y(1);y(2)], oeC0, dOm0, dt, P);
    try
        y = fmincon(fun24, xstar(2:3), [],[],[],[], [-0.2;-2*pi], [0.2;2*pi], nonl24, opts);
        x = [xfix1; y(:)];
    catch
        x = [xfix1; xstar(2); xstar(3)];
    end

    dv_ai = eor_cost_aiOmega(x, a, inc, da0, di0, V) * V;
    dv_e = eccentricity_corrected_cost(x, oeC0, oeT, dt, P, V, da0, di0);
    dvN = max(dv_ai, dv_e); % conservative, keeps eccentricity correction from reducing base cost
    if dvN < dvbest, dvbest = dvN; end
end

if ~isfinite(dvbest)
    dvbest = eor_cost_aiOmega(xstar, a, inc, da0, di0, V) * V;
end
dvopt = dvbest;
end

function f = eor_cost_aiOmega(x, a, inc, da0, di0, V)
da = a*x(1); di=x(2); dOm=x(3);
dv0 = sqrt((da/(2*a))^2 + di^2 + (dOm*sin(inc)/2)^2);
dvf = sqrt(((da-da0)/(2*a))^2 + (di-di0)^2 + (dOm*sin(inc)/2)^2);
f = dv0 + dvf; % normalized by V
end

function [c,ceq] = con_aiOmega(x, oeC0, dOm0, dt, P)
a=oeC0(1); inc=oeC0(3);
OmBase = rates_j2(oeC0,P); OmBase=OmBase(1);
oeShift = oeC0; oeShift(1)=a*(1+x(1)); oeShift(3)=inc+x(2);
OmShift = rates_j2(oeShift,P); OmShift=OmShift(1);
dOmDrift = (OmShift-OmBase)*dt;
ceq = wrapToPiLocal(x(3) + dOmDrift - dOm0)/(max(abs(dt*OmBase),1e-12));
c = [];
end

function xfix = phase_fixed_x1(xstar, oeC0, oeT, dt, P, N)
a=oeC0(1);
oeShift=oeC0; oeShift(1)=a*(1+xstar(1)); oeShift(3)=oeC0(3)+xstar(2);
ratesBase = rates_j2(oeC0,P);
ratesShift = rates_j2(oeShift,P);
oeCtf = propagate_oe_j2(oeC0, dt, P);
du0 = wrapToPiLocal((oeT(5)+oeT(6)) - (oeCtf(5)+oeCtf(6)));
du = wrapToPiLocal(du0 - ((ratesShift(2)-ratesBase(2)) + (ratesShift(3)-ratesBase(3)))*dt);
n = sqrt(P.mu/a^3);
daphase = -2*a*(du + 2*N*pi)/(3*n*dt);
xfix = xstar(1) + daphase/a;
end

function dv = eccentricity_corrected_cost(x, oeC0, oeT, dt, P, V, da0, di0)
a=oeC0(1); inc=oeC0(3);
da=a*x(1); di=x(2); dOm=x(3);
oeShift=oeC0; oeShift(1)=a+da; oeShift(3)=inc+di;
ratesShift = rates_j2(oeShift,P);
omega = oeC0(5) + ratesShift(2)*dt;
dex = oeT(2)*cos(oeT(5)) - oeC0(2)*cos(omega);
dey = oeT(2)*sin(oeT(5)) - oeC0(2)*sin(omega);

u0 = atan2(di, dOm*sin(inc)/2);
uf = atan2((di0-di), dOm*sin(inc)/2);
u0p = u0 + ratesShift(2)*dt;

base0 = (da/(2*a))^2 + di^2 + (dOm*sin(inc)/2)^2;
basef = ((da-da0)/(2*a))^2 + (di-di0)^2 + (dOm*sin(inc)/2)^2;

Aeq = [cos(u0p) cos(uf) sin(u0p) sin(uf);
       sin(u0p) sin(uf) -cos(u0p) -cos(uf)];
beq = [dex; dey];

fun = @(k) V*(sqrt(base0 + k(3)^2 + max(0,k(1))^2) + ...
              sqrt(basef + k(4)^2 + max(0,k(2))^2));
opts = optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',80);
try
    k = fmincon(fun, zeros(4,1), [],[], Aeq, beq, [], [], [], opts);
    dv = fun(k);
catch
    dv = V*(sqrt(base0)+sqrt(basef)) + V*norm([dex dey]);
end
end

function y=wrapToPiLocal(x)
y = mod(x+pi,2*pi)-pi;
end
