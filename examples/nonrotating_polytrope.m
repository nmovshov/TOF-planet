%% INTERIOR STRUCTURE OF A NON-ROTATING INDEX 1 POLYTROPE
% Example and test of the TOFPlanet class. We construct and converge a
% model of a non-rotating fluid planet with the pressure-density law
%
% $$P = K\rho^2$$
%
% with a polytropic constant $K$. This script demonstrates how to set up a
% TOFPlanet object with a specified barotrope and converge to a density
% structure in hydrostatic equilibrium with the given barotrope. The default
% starting density profile is that of a homogeneous sphere. This scripts also
% explores the error due to the radius grid resolution.
%
% Reminder: the density as a function of radius of a non-rotating, index-1
% polytropic planet in hydrostatic equilibrium is
% 
% $$\rho(r) = \rho_c \sin(ar)/(ar)$$
%
% where
% 
% $$a = \sqrt{2\pi{G}/K}$$
% 
% and
% 
% $$\rho_c = 3M/(4\pi{R^3})$$
%
% where $M$ is the planet's mass and $R = \pi/a$, which is independent of mass.
% 

%% Prepare workspace
clear
clc
close all
debug = true;
if debug
    u = setUnits;
else
    u = setFUnits;
end
G = u.gravity;

%% Set up some TOF Planets with arbitrary mass and radius
M = 317.8*u.earth_mass;
R = 71492*u.km;

%#ok<*SAGROW>
N = [128,256,512];
for k=1:length(N)
    tof(k) = TOFPlanet('debug', debug);
    tof(k).name = ['TOF',int2str(N(k))];
    tof(k).mass = M;
    tof(k).radius = R;
    tof(k).si = R*linspace(1, 1/N(k), N(k))';
    tof(k).rhoi = ones(N(k),1)*M/(4*pi/3*R^3);
    tof(k).mrot = 0;
    tof(k).P0 = 0.0*u.bar;
end

%% Construct a polytrope of index 1 to represent the planet's eos
% A non-rotating, index-1 polytrope is completely defined by K. The radius is
% independent of mass. The density structure normalized to the central density
% is also independent of mass, and the absolute value of density is determined
% with the known average density.
n = 1; % polytrope index
K = 2*G/pi*R^2; % polytrope constant
eos = barotropes.Polytrope(K, n);
eos.name = sprintf('$P\\propto\\rho^2$');
[tof.eos] = deal(eos);

%% Relax to desired barotrope
for k=1:length(N)
    tof(k).opts.verbosity = 1;
    tof(k).opts.dJtol = 1e-6;
    tof(k).opts.drhotol = 1e-6;
    tof(k).opts.MaxIterBar = 20; % default 40
    tof(k).relax_to_barotrope;
end

%% Compare computed and analytic density structure
% prepare
set(groot, 'defaultTextInterpreter', 'latex')
set(groot, 'defaultLegendInterpreter', 'latex')
set(groot, 'defaultAxesBox', 'on')

% calculate
a = sqrt(2*pi*G/K);
R = pi/a;
rho_av = 3*M/(4*pi*R^3);
rho_c = (pi^2/3)*rho_av;
r = linspace(0,1)*R;
rho_exact = rho_c*sin(a*r)./(a*r);
rho_exact(1) = rho_c;

% plot
ah = axes; hold(ah, 'on');
l1 = plot(r/R, rho_exact/rho_c, 'k--');
l1.DisplayName = '$\sin(ar)/(ar)$';
for k=1:length(N)
    l2(k) = plot(tof(k).si/R, tof(k).rhoi/rho_c, '-');
    l2(k).DisplayName = ['$\rho_i/\rho_c$ ',tof(k).name];
end

% annotate
xlabel('$r/R$')
ylabel('$\rho/\rho_c$')

lh = legend(ah, 'show');
lh.FontSize = 12;

% errors
n = length(N);
P_center = interp1(double(tof(n).si), double(tof(n).Pi), 0, 'pchip')*u.Pa;
P_err = (P_center - tof(n).eos.pressure(rho_c))/tof(n).eos.pressure(rho_c);
s_tit = sprintf(['$N=%g$; ',...
    '$\\beta=%g$; $\\Delta P_c=%g\\%%$'],...
    N(n), double(tof(n).betam), double(P_err)*100);
title(s_tit)

% polytrope definitions
s_poly{1} = '$a=\sqrt{2\pi{}G/K}$';
s_poly{2} = '$R=\pi/a$';
s_poly{3} = '$\rho_c = \frac{\pi{M}}{4R^3}$';
th = annotation('textbox');
th.Interpreter = 'latex';
th.String = s_poly;
th.FontSize = 12;
th.Position = [lh.Position(1), lh.Position(2)-th.Position(4),...
    lh.Position(3), th.Position(4)];

% Also plot converged and input barotrope
tof(end).plot_barotrope('showinput',true,'showscaledinput',true);
title(s_tit);
