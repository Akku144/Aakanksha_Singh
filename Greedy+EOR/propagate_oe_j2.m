function oe2 = propagate_oe_j2(oe, dt, P)
rates = rates_j2(oe, P);
oe2 = oe;
oe2(4) = wrapTo2PiLocal(oe(4) + rates(1)*dt);
oe2(5) = wrapTo2PiLocal(oe(5) + rates(2)*dt);
oe2(6) = wrapTo2PiLocal(oe(6) + rates(3)*dt);
end
function y = wrapTo2PiLocal(x)
y = mod(x,2*pi);
end
