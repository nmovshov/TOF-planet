%% test convergence of NMoI
clear
clc

zstrat = @zvecs.exponential;
name = char(zstrat);
name = name(7:end);

x = 2.^(7:14);
for k=1:length(x)
    tof = tofmodels.triple_polytrope(x(k), [1.14e5, 1, 4.6e4, 0.9, 30, 0.67, 0.85, 0.1], zstrat);
    tof.mass = 2e27;
    tof.radius = 7e7;
    tof.mrot = 0;
    tof.opts.verbosity = 2;
    tof.opts.MaxIterBar = 60;
    
    tof.relax_to_barotrope;
    L(k) = tof.NMoI;
    beta(k) = tof.betam;
    clear tof
end
clear k
save
