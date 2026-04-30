function P = constants_eor()
P.mu = 398600.4418;       % km^3/s^2
P.Re = 6378.173;          % km, as used in Ye et al.
P.J2 = 1.08262668e-3;
P.amin = P.Re + 150;      % paper only says limited by Earth radius; choose safe LEO lower bound
P.day = 86400;
end
