%% Debugging TOFPlanet
clear
clc
close all
si = setFUnits;

N = 128;
tof = TOFPlanet(N);
tof.G = si.gravity;
tof.M = 2e27*si.kg;
tof.Rm = 7e7*si.m;
tof.rhoi = ones(N,1)*tof.M/(4*pi/3*tof.Rm^3);

poly_n = 1;
poly_K = 2e5*si.Pa/(si.kg/si.m^3)^(1 + 1/poly_n);
eos1 = barotropes.Polytrope(poly_K, poly_n);
tof.eos = eos1;

tof.opts.verbosity = 1;
%eos0 = barotropes.ConstDensity(eos1.density(1*si.bar));
%tof.eos = [eos0; repmat(eos1, N-1, 1)];
