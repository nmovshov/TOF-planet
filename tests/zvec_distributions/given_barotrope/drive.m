%% test convergence 
clear
clc

zstrat = @zvecs.cosine;
name = char(zstrat);
name = name(7:end);

x = 2.^(7:14);
for k=1:length(x)
    tof = tofmodels.triple_polytrope(x(k), [1.14e5, 1, 4.6e4, 0.9, 30, 0.67, 0.85, 0.1], zstrat);
    tof.mass = 2e27;
    tof.radius = 7e7;
    tof.mrot = 0.1;
    tof.opts.verbosity = 2;
    tof.opts.MaxIterBar = 60;
    tof.opts.dJtol = 1e-8;
    
    tof.relax_to_barotrope;
    L(k) = tof.NMoI;
    beta(k) = tof.betam;
    y2(k) = tof.J2;
    y4(k) = tof.J4;
    y6(k) = tof.J6;
    clear tof
end
clear k ans
save
