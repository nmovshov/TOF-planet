%% Test tof7 barotrope convergence with varying xlevels
% This script runs the n=1 polytrope with large N and varying xlevels to test
% convergence and runtime performance of skip-n-spline with tof7.
%
% Convergence is gauged against the full layer count case, but just for fun we 
% aim for the numbers of Wisdom and Hubbard (2016) Table 3.

%% Prepare workspace
clear
clc
close all
u = setFUnits;

%% Construct a polytrope of index 1, aiming for exact replicaiton of WH16
G = 6.6738480e-11; % Hubbard to Guillot to me personal communcation
GM = 1.266865361e17; % WH16
M = GM/G;
Re = 71492*u.km; % (to match K use K=2*G/pi*R^2 instead, but it doesn't matter)
qrot = 0.089195487; % WH16
wrot = sqrt(qrot*GM/Re^3);
Prot = 2*pi/wrot;
aos = 1.022875431133185; % WH16 Table 3 Re/R
K = 2.003565e5; % Hubbard to Guillot personal communication (no effect on Js)
n = 1;
eos = barotropes.Polytrope(K, n);
eos.name = '$P\propto\rho^2$';

%% Set up TOFPlanet(s)
N = 2^16;
nx = 2.^(4:16);

for k=1:length(nx)
    tof = TOFPlanet('toforder',7,'xlevels',nx(k));
    tof.name = [int2str(nx(k)),'-xlevels'];
    tof.G = G; % undocumented TOFPlanet property
    tof.mass = M;
    tof.radius = Re;
    tof.period = Prot; % trying to match WH16 qrot
    tof.si = Re*linspace(1, 1/N, N)'; % will be renormalized
    tof.rhoi = ones(N,1)*M/(4*pi/3*Re^3); % will be renormalized
    tof.P0 = 0*u.bar; % added to surface pressure
    tof.eos = eos;
    TOFS(k) = tof;
end

%% Relax to desired barotrope (fast for tof4, slow for tof7)
textprogressbar('Running ToF7 models...')
for k=1:length(nx)
    tof = TOFS(k);
    tof.opts.verbosity = 0;
    tof.opts.drhotol = 1e-6;
    tof.opts.dJtol = 1e-10;
    tof.opts.MaxIterBar = 60;
    tof.opts.MaxIterHE = 60;
    rtime(k) = tof.relax_to_barotrope();
    textprogressbar(k/length(N)*100)
end
textprogressbar(' done.')

%% Construct the benchmarking table
% The variables to compare are [Re/R, J2, J4, ..., J14]

% Wisdom and Hubbard (2016) Table 3
CLC = [nan, 1.022875431133185, 1.398851089834637e-2, -5.318281001092471e-4,...
                          3.011832290533577e-5, -2.132115710726158e-6,...
                          1.740671195871128e-7, -1.568219505602588e-8,...
                          1.518099230068580e-9];

% With TOFPlanet
for k=1:length(nx)
    tof = TOFS(k);
    MTOF7(k,:) = [nx(k), rtime(k), tof.a0/tof.s0, tof.Js(2:end)];
end

cols = {'nx', 'runtime', 'Re/R', 'J2', 'J4', 'J6', 'J8', 'J10', 'J12', 'J14'};
rows = {TOFS.name};
A = MTOF7;
E = (A(1:end-1,3:6) - A(end,3:6))./A(end,3:6);
T_vals = array2table(A, 'VariableNames', cols, 'RowNames', rows);
T_errs = array2table(E, 'VariableNames', cols(3:6), 'RowNames', rows(1:end-1));

%% Output
format shorte
display(T_errs)
format
save n1ofnx.mat
