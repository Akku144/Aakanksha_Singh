clear; clc;

P = constants_eor();
[targets, ssc] = data_40targets();

cfg.Ns = 2;
cfg.Nt = 40;
cfg.NmaxPhase = 1;

seq = cell(2,1);

seq{1} = [30 6 16 28 33 24 20 19 2 11 26 29 18 1 12 14 39 25 17 31 21 7 22 35 23 37 13 40 34 9 5 32];
seq{2} = [36 3 4 27 15 10 8 38];

tic;
J = sequence_cost_direct(ssc, targets, seq, P, cfg);
compTime = toc;

fprintf('Paper sequence fuel from our EOR = %.6f km/s\n', J);
fprintf('Computational time = %.3f s\n', compTime);

function J = sequence_cost_direct(ssc0, targets, seq, P, cfg)

J = 0;
ssc = ssc0;
tnow = zeros(size(ssc0,1),1);

for k = 1:numel(seq)
    for j = seq{k}
        dt = transfer_time(ssc(k,:), targets(j,:), P);
        dv = eor_estimate(ssc(k,:), targets(j,:), dt, P, cfg.NmaxPhase);

        J = J + dv;
        tnow(k) = tnow(k) + dt;
        ssc(k,:) = propagate_oe_j2(targets(j,:), tnow(k), P);
    end
end

end