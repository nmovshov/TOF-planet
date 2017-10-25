%% GRAVITY COEFFICIENTS OF ROTATING INDEX-1 POLYTROPE
% Example and test of the TofPlanet class. We construct and converge a model of a
% rotating fluid planet with the pressure-density law
%
% $$P = K\rho^2$$
%
% with a polytropic constant $K$. This script demonstrates how to set up a
% TOFPlanet object with a specified barotrope and converge to a density structure
% in hydrostatic equilibrium with the given barotrope. The default starting
% density profile is that of a homogeneous sphere.
%
% The comparison model is the 3rd-order Zharkov and Trubitsyn theory. I take the
% numbers from Table 5 of Hubbard (2013).

%% Prepare workspace
clear
clc
close all
debug = false;
if debug
    u = setUnits;
else
    u = setFUnits;
end
G = u.gravity;

%% Set up a TOF Planet with arbitrary mass and radius
% I'm using numbers for Jupiter just for kicks, no effect on Js of course.
M = 317.8*u.earth_mass;
%R = 6.9917979e7*u.m;
R = 71492*u.km;

N = 1024;
tof = TOFPlanet(N,'debug',debug);
tof.name = [int2str(N),'-point TOF'];
tof.mass = M;
tof.radius = R;
tof.si = R*linspace(1, 1/N, N)'; % will be renormalized
tof.rhoi = ones(N,1)*M/(4*pi/3*R^3); % will be renormalized
tof.mrot = 0.083432862292087;
tof.P0 = 0*u.bar; % added to surface pressure

%% Construct a polytrope of index 1 to represent the planet's eos
n = 1;
K = 2*G/pi*R^2; % ...matches radius just for show, K has no effect on the Js
eos = barotropes.Polytrope(K, n);
eos.name = '$P\propto\rho^2$';
tof.eos = eos;

%% To (barely) speed up convergence start with an approximate density structure
a = sqrt(2*pi*G/K);
r = pi/a;
rho_av = 3*M/(4*pi*r^3);
rho_c = (pi^2/3)*rho_av;
x = (tof.si(1:end-1) + tof.si(2:end))/2;
x(end+1) = tof.si(end)/2;
tof.rhoi = rho_c*sin(a*x)./(a*x);

%% Relax to desired barotrope
tof.opts.drhotol = 1e-5;
tof.opts.dJtol = 1e-6;
tof.opts.MaxIterBar = 30;
tof.relax_to_barotrope;

%% Compare computed and analytic density structure
q = tof.qrot;
% Zharkov & Trubistyn (1978) eq. 34.12
ZT3 = [q;...
    (0.173273*q - 0.197027*q^2 + 0.15*q^3)*1e2;...
    (-0.081092*q^2 + 0.15*q^3)*-1e4;...
    (0.056329*q^3)*1e5;...
    nan; nan; nan];
% Hubbard (2013) Table 5
H13_256 = [q; 1.3991574; 5.3182810; 3.0118323; 2.1321157; 1.7406710; 1.5682179];
H13_512 = [q; 1.3989253; 5.3187997; 3.0122356; 2.1324628; 1.7409925; 1.5685327];

% TOFPlanet
TOF = [q; tof.J2*1e2; -tof.J4*1e4; tof.J6*1e5; -tof.J8*1e6; nan; nan];

% Make it a table
T = table(ZT3, H13_256, H13_512, TOF);
T.Properties.RowNames = {'q','J2x10^2','-J4x10^4','J6x10^5','-J8x10^6',...
    'J10x10^7','-J12x10^8'};

% Display
format long
format compact
disp(T)
format
format compact
J2_err = (tof.J2*1e2 - H13_512(2))/(tof.J2*1e2)
J4_err = (-tof.J4*1e4 - H13_512(3))/(-tof.J4*1e4)
J6_err = (1e5*tof.J6 - H13_512(4))/(tof.J6*1e5)
format
try
    %tof.plot_equipotential_surfaces;
    tof.plot_barotrope('showinput',true,'showscaledinput',true);
catch
end
