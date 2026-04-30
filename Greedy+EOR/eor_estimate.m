function dvopt = proxy_cost(oeC0, oeT0, dt, P, Nmax)
% Clean EOR estimator structure based on Huang Eq. (14)-(34).

if nargin < 5, Nmax = 5; end
if dt <= 0 || ~isfinite(dt), dvopt = inf; return; end

oeC = propagate_oe_j2(oeC0, dt, P);
oeT = propagate_oe_j2(oeT0, dt, P);

a = oeC0(1); inc = oeC0(3); V = sqrt(P.mu/a);
rates0 = rates_j2(oeC0, P); Om0dot = rates0(1);

dOm0 = wrapToPiLocal(oeT(4) - oeC(4));
da0 = oeT(1) - oeC(1);
di0 = oeT(3) - oeC(3);

den = 49/2 + tan(inc)^2/2 + 2/(sin(inc)^2*(dt*Om0dot)^2);
lam = (dOm0/(dt*Om0dot)) / den;
xseed = [7*lam; 0.5*tan(inc)*lam; -2*lam/(sin(inc)^2*dt*Om0dot)];

opts_fmin = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
    'MaxIterations',100,'OptimalityTolerance',1e-10,'ConstraintTolerance',1e-10);

fun19 = @(x) eor_cost_aiOmega(x, a, inc, da0, di0);
nonl19 = @(x) con_aiOmega(x, oeC0, dOm0, dt, P);
lb = [(P.amin-a)/a; -0.2; -2*pi];
ub = [0.2; 0.2; 2*pi];

try
    xstar = fmincon(fun19, xseed, [], [], [], [], lb, ub, nonl19, opts_fmin);
catch
    xstar = xseed;
end

dvbest = inf;

for N = -Nmax:Nmax
    xfix1 = phase_fixed_x1(xstar, oeC0, oeT, dt, P, N);

    if xfix1 < (P.amin-a)/a || xfix1 > 0.2
        continue;
    end

    fun24 = @(y) eor_cost_aiOmega([xfix1; y(1); y(2)], a, inc, da0, di0);
    nonl24 = @(y) con_aiOmega([xfix1; y(1); y(2)], oeC0, dOm0, dt, P);

    try
        y = fmincon(fun24, xstar(2:3), [], [], [], [], [-0.2; -2*pi], [0.2; 2*pi], nonl24, opts_fmin);
        x = [xfix1; y(:)];
    catch
        x = [xfix1; xstar(2); xstar(3)];
    end

    dv_e = eccentricity_corrected_cost(x, oeC0, oeT0, dt, P, V, da0, di0);

    if dv_e < dvbest
        dvbest = dv_e;
    end
end

if ~isfinite(dvbest)
    dvbest = V * eor_cost_aiOmega(xstar, a, inc, da0, di0);
end

dvopt = dvbest;

%scale = 0.5272;   % calibration to Ye Table 11 paper sequence
%dvopt = scale * dvbest;
end

function f = eor_cost_aiOmega(x, a, inc, da0, di0)
da = a*x(1); di = x(2); dOm = x(3);
dv0 = sqrt((da/(2*a))^2 + di^2 + (dOm*sin(inc)/2)^2);
dvf = sqrt(((da-da0)/(2*a))^2 + (di-di0)^2 + (dOm*sin(inc)/2)^2);
f = dv0 + dvf;
end

function [c, ceq] = con_aiOmega(x, oeC0, dOm0, dt, P)
a = oeC0(1); inc = oeC0(3);
OmBase = rates_j2(oeC0, P); OmBase = OmBase(1);
oeShift = oeC0;
oeShift(1) = a*(1 + x(1));
oeShift(3) = inc + x(2);
OmShift = rates_j2(oeShift, P); OmShift = OmShift(1);
dOmDrift = (OmShift - OmBase)*dt;
ceq = (x(3) + dOmDrift - dOm0)/(max(abs(dt*OmBase),1e-12));
c = [];
end

function xfix = phase_fixed_x1(xstar, oeC0, oeTtf, dt, P, N)
a = oeC0(1);
oeShift = oeC0;
oeShift(1) = a*(1 + xstar(1));
oeShift(3) = oeC0(3) + xstar(2);
ratesBase = rates_j2(oeC0, P);
ratesShift = rates_j2(oeShift, P);
oeCtf = propagate_oe_j2(oeC0, dt, P);
du0 = wrapToPiLocal((oeTtf(5) + oeTtf(6)) - (oeCtf(5) + oeCtf(6)));
du = wrapToPiLocal(du0 - ((ratesShift(2)-ratesBase(2)) + (ratesShift(3)-ratesBase(3)))*dt);
n = sqrt(P.mu/a^3);
daphase = -2*a*(du + 2*N*pi)/(3*n*dt);
xfix = xstar(1) + daphase/a;
end

function dv = eccentricity_corrected_cost(x, oeC0, oeT0, dt, P, V, da0, di0)
a = oeC0(1); i = oeC0(3);
da = a*x(1); di = x(2); dOmega = x(3);

oeShift = oeC0;
oeShift(1) = a + da;
oeShift(3) = i + di;
ratesShift = rates_j2(oeShift, P);

omega = oeC0(5) + ratesShift(2)*dt;
oeTtf = propagate_oe_j2(oeT0, dt, P);

dex = oeTtf(2)*cos(oeTtf(5)) - oeC0(2)*cos(omega);
dey = oeTtf(2)*sin(oeTtf(5)) - oeC0(2)*sin(omega);

u0 = atan2(di, dOmega*sin(i)/2);
uf = atan2((di0-di), dOmega*sin(i)/2);
u0p = u0 + ratesShift(2)*dt;

Aeq = [cos(u0p), cos(uf),  sin(u0p),  sin(uf);
       sin(u0p), sin(uf), -cos(u0p), -cos(uf)];
beq = [dex; dey];

opts_fsolve = optimoptions('fsolve','Display','off','FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12,'MaxIterations',200);

bestCase = inf;

for caseID = 1:4
    F = @(z) kkt_system(z, caseID, Aeq, beq, da, da0, di, di0, dOmega, i, V, a);
    k0 = Aeq \ beq;
    z0 = [k0; 0; 0];

    try
        zsol = fsolve(F, z0, opts_fsolve);
        k = zsol(1:4);

        if case_valid(k, caseID, da, da0, a)
            val = eccentricity_case_cost(k, caseID, da, da0, di, di0, dOmega, i, V, a);
            bestCase = min(bestCase, val);
        end
    catch
    end
end

if ~isfinite(bestCase)
    opts_fmin = optimoptions('fmincon','Display','off','Algorithm','sqp');
    for caseID = 1:4
        fun = @(k) eccentricity_case_cost(k, caseID, da, da0, di, di0, dOmega, i, V, a);
        try
            ksol = fmincon(fun, Aeq\beq, [], [], Aeq, beq, [], [], [], opts_fmin);
            bestCase = min(bestCase, fun(ksol));
        catch
        end
    end
end

dv = bestCase;
end

function tf = case_valid(k, caseID, da, da0, a)
k1 = k(1); k2 = k(2);
b1 = da/a; b2 = (da-da0)/a;
tol = 1e-8;
switch caseID
    case 1
        tf = (k1 <= b1 + tol) && (k2 <= b2 + tol);
    case 2
        tf = (k1 <= b1 + tol) && (k2 >= b2 - tol);
    case 3
        tf = (k1 >= b1 - tol) && (k2 <= b2 + tol);
    case 4
        tf = (k1 >= b1 - tol) && (k2 >= b2 - tol);
    otherwise
        tf = false;
end
end

function cost = eccentricity_case_cost(k, caseID, da, da0, di, di0, dOmega, i, V, a)
k1 = k(1); k2 = k(2); k3 = k(3); k4 = k(4);
A0 = di^2 + (dOmega*sin(i)/2)^2;
Af = (di-di0)^2 + (dOmega*sin(i)/2)^2;

switch caseID
    case 1
        dv0 = V*sqrt((da/(2*a))^2 + A0 + k3^2);
        dvf = V*sqrt(((da-da0)/(2*a))^2 + Af + k4^2);
    case 2
        dv0 = V*sqrt((da/(2*a))^2 + A0 + k3^2);
        dvf = V*sqrt((k2/2)^2 + Af + k4^2);
    case 3
        dv0 = V*sqrt((k1/2)^2 + A0 + k3^2);
        dvf = V*sqrt(((da-da0)/(2*a))^2 + Af + k4^2);
    case 4
        dv0 = V*sqrt((k1/2)^2 + A0 + k3^2);
        dvf = V*sqrt((k2/2)^2 + Af + k4^2);
end
cost = dv0 + dvf;
end

function F = kkt_system(z, caseID, Aeq, beq, da, da0, di, di0, dOmega, i, V, a)
k = z(1:4);
lambda = z(5); chi = z(6);
grad = numerical_grad(@(kk) eccentricity_case_cost(kk, caseID, da, da0, di, di0, dOmega, i, V, a), k);
g = Aeq*k - beq;
F = zeros(6,1);
F(1:4) = grad + lambda*Aeq(1,:).' + chi*Aeq(2,:).';
F(5:6) = g;
end

function g = numerical_grad(fun, x)
h = 1e-7;
g = zeros(size(x));
for j = 1:numel(x)
    xp = x; xm = x;
    xp(j) = xp(j) + h;
    xm(j) = xm(j) - h;
    g(j) = (fun(xp) - fun(xm))/(2*h);
end
end
