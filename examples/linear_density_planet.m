%% GRAVITY COEFFICIENTS OF LINEAR DENSITY PLANET
% Example and test of the TOFPlanet class. We construct and converge a
% model of a rotating fluid planet with a density profile that is linear in the
% mean radius of a level surface.
%
% The comparison is with a 5th-order Zharkov and Trubitsyn model. I take the
% numbers from Table 1 of Hubbard (2013).

%% Prepare workspace
clear
clc
close all

%% Set up a TOF object and give it a density profile linear in mean radii
N = 2048;
tof = TOFPlanet();
tof.si = linspace(1, 1/N, N)';
dl = 1/(N-1);
tof.rhoi = [0:dl:1]';
tof.mrot = 0.0830; % Hubbard (2013) Table 1

%% Relax to hydrostatic equilibrium
tof.opts.verbosity = 2;
tof.opts.dJtol = 1e-12;
tof.relax_to_HE;

%% Compare computed and analytic density structure
q = tof.qrot;
m = tof.mrot;

% Zharkov & Trubistyn (1978) Table 3.1
ZT5 = [nan; 0.0830; 1.4798; 5.929; 3.497; 2.52; 2.4; nan; nan];
% Hubbard (2013) Table 1
H13_128 = [0.088822426; 0.082999915; 1.4798138; 5.9269129; 3.4935680; 2.5493209;...
    2.1308951; 1.9564143; 1.9237724];

% TOFPlanet
TOF = [q; m; tof.J2*1e2; -tof.J4*1e4; tof.J6*1e5; -tof.J8*1e6; nan; nan; nan];

% Make it a table
T = table(ZT5, H13_128, TOF);
T.Properties.RowNames = {'q','m','J2x10^2','-J4x10^4','J6x10^5','-J8x10^6',...
    'J10x10^7','-J12x10^8','J14x10^9'};

% Display
format long
format compact
fprintf('\n')
disp(T)
format
