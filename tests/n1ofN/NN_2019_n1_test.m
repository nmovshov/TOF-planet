%% Test tof4 and tof7 convergence with N on the n=1 polytrope
% This script replicates Nadine's tof benchmarks against the CLC solution of
% WH16. The n=1 polytrope is run with tof<n> with increasing number of level
% radii. I won't be using skip-n-spline since Nadine didn't and I want a fair
% comparison.
%
% The numbers to benchmark against are apparently those of the CLC method
% (Wisdom, 1996, unpublished) which is lost to time. I take the numbers from
% Wisdom and Hubbard (2016) Table 3. For consistency we try to replicate the
% Wisdom and Hubbard model exactly, matching their q value exactly by using
% also their G and K values (Guillot, double hearsay).

%% Prepare workspace
clear
clc
close all

%% Choose layer numbers and distribution strategies to investigate
N = [128, 256];
nx = -1;
zstrat = @(n)zvecs.topheavy(n);

%% Construct a polytrope of index 1, aiming for exact replicaiton of WH16
G = 6.6738480e-11; % Hubbard to Guillot to me personal communcation
GM = 1.266865361e17; % WH16
M = GM/G;
Re = 71492*1e3; % (to match K use K=2*G/pi*R^2 instead, but it doesn't matter)
qrot = 0.089195487; % WH16
wrot = sqrt(qrot*GM/Re^3);
Prot = 2*pi/wrot;
aos = 1.022875431133185; % WH16 Table 3 Re/R
K = 2.003565e5; % Hubbard to Guillot personal communication (no effect on Js)
n = 1;
eos = barotropes.Polytrope(K, n);
eos.name = '$P\propto\rho^2$';

%% Set up the TOFPlanet(s)
for k=1:length(N)
    tofour = TOFPlanet('toforder',4);
    tofour.name = sprintf('%d-%s TOF4',...
        N(k),string(char(zstrat)).extractBetween('.','('));
    tofour.G = G; % undocumented TOFPlanet property
    tofour.mass = M;
    tofour.radius = Re;
    tofour.period = Prot; % trying to match WH16 qrot
    tofour.si = Re*zstrat(N(k)); % will be renormalized
    tofour.rhoi = ones(N(k),1)*M/(4*pi/3*Re^3); % will be renormalized
    tofour.P0 = 0.0; % added to surface pressure
    tofour.eos = eos;
    FOURS(k) = tofour;
    
    tofsev = TOFPlanet('toforder',7);
    tofsev.name = sprintf('%d-%s TOF7',...
        N(k),string(char(zstrat)).extractBetween('.','('));
    tofsev.G = G; % undocumented TOFPlanet property
    tofsev.mass = M;
    tofsev.radius = Re;
    tofsev.period = Prot; % trying to match WH16 qrot
    tofsev.si = Re*zstrat(N(k)); % will be renormalized
    tofsev.rhoi = ones(N(k),1)*M/(4*pi/3*Re^3); % will be renormalized
    tofsev.P0 = 0.0; % added to surface pressure
    tofsev.eos = eos;
    SEVENS(k) = tofsev;
end

%% Relax to desired barotrope (fast for tof4, slow for tof7)
textprogressbar('Running ToF4 models...') 
t = tic;
for k=1:length(N)
    tofour = FOURS(k);
    tofour.opts.xlevels = nx;
    tofour.opts.verbosity = 0;
    tofour.opts.drhotol = 1e-6;
    tofour.opts.dJtol = 1e-10;
    tofour.opts.MaxIterBar = 60;
    tofour.opts.MaxIterHE = 60;
    rt4(k) = tofour.relax_to_barotrope();
    textprogressbar(k/length(N)*100)
end
textprogressbar(sprintf(' done. (%s)',seconds2human(toc(t))))

textprogressbar('Running ToF7 models...')
t = tic;
for k=1:length(N)
    tofsev = SEVENS(k);
    tofsev.opts.xlevels = nx;
    tofsev.opts.verbosity = 0;
    tofsev.opts.drhotol = 1e-6;
    tofsev.opts.dJtol = 1e-10;
    tofsev.opts.MaxIterBar = 60;
    tofsev.opts.MaxIterHE = 60;
    rt7(k) = tofsev.relax_to_barotrope();
    textprogressbar(k/length(N)*100)
end
textprogressbar(sprintf(' done. (%s)',seconds2human(toc(t))))

%% Construct the benchmarking table
% The variables to compare are [Re/R, J2, J4, ..., J14]

% Wisdom and Hubbard (2016) Table 3
CLC = [nan, nan, 1.022875431133185,...
        1.398851089834637e-2, -5.318281001092471e-4, 3.011832290533577e-5,...
       -2.132115710726158e-6, 1.740671195871128e-7, -1.568219505602588e-8,...
        1.518099230068580e-9];

% With TOFPlanet
for k=1:length(N)
    tofour = FOURS(k);
    tofsev = SEVENS(k);
    MTOF4(k,:) = [N(k),nx,tofour.a0/tofour.s0,tofour.Js(2:end),nan,nan,nan];
    MTOF7(k,:) = [N(k),nx,tofsev.a0/tofsev.s0,tofsev.Js(2:end)];
end

cols = {'N', 'nx', 'Re/R', 'J2', 'J4', 'J6', 'J8', 'J10', 'J12', 'J14'};
rows = {'CLC', FOURS.name, SEVENS.name};
A = [CLC; MTOF4; MTOF7];
E = (A(2:end,:) - A(1,:))./A(1,:);
E(:,1:2) = A(2:end,1:2);
T_vals = array2table(A, 'VariableNames', cols, 'RowNames', rows);
T_errs = array2table(E, 'VariableNames', cols, 'RowNames', rows(2:end));

%% Output and save
format shorte
display(T_errs(:,3:4))
format
clear SEVENS FOURS tofour tofsev k n u
save(sprintf('%f.mat',now()), 'T_vals','T_errs')
