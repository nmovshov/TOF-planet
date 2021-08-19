%% An example of a three-layer planet, with density jumps
% As a simple example of a planet that exhibits desnity jumps I'll use a
% parameterization based on 8th-degree polynomial with arctan smoothed step
% functions.

%% Prepare workspace
clear
clc
close all

%% Load parameter values for the ppwd representation and generate rho(s)
load ppwd_uranus.mat % this one approximates model U2 of Nettelmann et al. 2013
N = 4096;
[zvec,dvec] = ppwd_profile(N,x,obs.rho0); % x and obs loaded from mat file

%% Set up a TOF object and give it the density profile with jumps
tof = TOFPlanet();
tof.name = 'PPWD-N13-U2';
tof.si = obs.s0*zvec;
tof.rhoi = dvec;
tof.set_observables(obs); % also sets tof.period.

%% Relax to hydrostatic equilibrium
tof.opts.xlevels = 256;
tof.relax_to_HE();
tof.fix_radius();

%% Examine planet
tof.plot_rho_of_r();
tof.report_card(obs)
