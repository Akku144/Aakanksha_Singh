function r = rates_j2(oe, P)
% returns [RAANdot, omegadot, Mdot], rad/s, near-circular J2 secular model
a=oe(1); e=oe(2); inc=oe(3);
p = a*(1-e^2);
n = sqrt(P.mu/a^3);
Omdot = -1.5*P.J2*(P.Re/p)^2*n*cos(inc);
wdot  =  1.5*P.J2*(P.Re/p)^2*n*(2 - 2.5*sin(inc)^2);
Mdot  =  n + 1.5*P.J2*(P.Re/p)^2*n*(1 - 1.5*sin(inc)^2)*sqrt(1-e^2);
r=[Omdot wdot Mdot];
end
