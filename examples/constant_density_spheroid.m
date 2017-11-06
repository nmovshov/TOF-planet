%% MACLAURIN ELLIPSOID WITH TOF
% The Maclaurin spheroid is a closed analytic solution for the shape of a
% rotating *constant density* self-gravitating fluid. This script compares the
% numerical ToF solution with the expected analytic solution.

%% Maclaurin's solution
% The equilibrium shape of a constant density rotating fluid can be shown
% (although not by me) to be an ellipsoid of revolution such that the radius $r$
% follows
% 
% $$r^2(\mu) = \frac{a^2}{1 + l^2\mu^2}$$
% 
% where $a$ is the equatorial radius, $\mu$ is the cosine of the angle from the
% rotation axis (the colatitude), $b$ is the polar radius and
%
% $$l^2 = \frac{a^2}{b^2} - 1.$$
%
% The ellipticity parameter $l$ is related to a dimensionless rotation
% parameter $m=\omega^2s^3/GM$ by the transcendental equation
%
% $$ m = \frac{3}{2l^3}[(3 + l^2)\arctan{l} - 3l]. $$
% 
% The rotation parameter $m$ is given in terms of the ellipsoid's *mean* radius
% $s$. For an oblate ellipsoid the radii are related by $s^3=ba^2$.

%% Prepare workspace
clear
clc
close all

%% Construct an exact normalized (a=1) Maclaurin ellipsoid
m = 0.1; % small rotation parameter
lfun = @(x)(3./(2*x.^3)).*((3 + x.^2).*atan(x) - 3*x) - m;
el = fzero(lfun,0.5); % ellipse parameter
b_exact = sqrt(1/(1 + el^2)); % polar radius
s3 = b_exact; % *mean* radius, s^3=b*a^2 (but a=1) we will use this later
xi_exact = @(mu)1./sqrt((1 + (el^2).*(mu.^2)));

%% Take a quick look for sanity check
try % requires R2016a or later
    theta = linspace(0,2*pi);
    polax = polarplot(theta, xi_exact(cos(theta)));
    polax.Parent.ThetaZeroLocation = 'top';
    polax.Parent.ThetaDir = 'clockwise';
    polax.Parent.ThetaAxisUnits = 'rad';
    hold
catch
end

%% Now set up a TOFPlanet to mimic constant density case
N = 12; % can be any number though!
tof = TOFPlanet();
tof.si = linspace(1,1/N,N)';
tof.rhoi = ones(N,1);
tof.mrot = m;
tof.opts.verbosity = 2;
tof.relax_to_HE();

%% Take a quick look for sanity check
try % requires R2016a or later
    theta = linspace(0,2*pi);
    polax = polarplot(theta, tof.level_surface(1, theta)/tof.level_surface(1,pi/2));
    polax.Parent.ThetaZeroLocation = 'top';
    polax.Parent.ThetaDir = 'clockwise';
    polax.Parent.ThetaAxisUnits = 'rad';
catch
end

%% Compare numerical and analytic solutions
% Compare the level surface radii
teta = 0:0.01:pi;
mu = cos(teta);
xi_tof = @(t)tof.level_surface(1,t);
dxi = xi_tof(teta)/xi_tof(pi/2) - xi_exact(mu);
figure
lh = semilogy(mu, abs(dxi));
ah = lh.Parent;
ah.XLabel.Interpreter = 'latex';
ah.YLabel.Interpreter = 'latex';
ah.Title.Interpreter = 'latex';
ah.XLabel.String = '$\mu = \cos(\theta)$';
ah.YLabel.String = '$d\xi$';
ah.Title.String = '$d\xi = \xi(\mu) - 1/\sqrt{1 + l^2\mu^2}$';
xi_err = max(abs(dxi)) %#ok<NOPTS>
%b_err = b_cms - b_exact %#ok<NOPTS>

% Compare the J values
n = 0:2:8;
J_exact = (-1).^(1 + n/2).*(3./((n + 1).*(n + 3))).*(el^2/(1 + el^2)).^(n/2);
J_tof = tof.Js;
dJ = J_tof - J_exact;
subplot(2,1,1,ah);
subplot(2,1,2);
lh = stem(n, abs(dJ));
ah = lh.Parent;
ah.YScale = 'log';
ah.XTick = n(1:2:end);
ah.XMinorTick = 'on';
ah.XLabel.Interpreter = 'latex';
ah.YLabel.Interpreter = 'latex';
ah.Title.Interpreter = 'latex';
ah.XLabel.String = '$n = 0,2,4,\ldots$';
ah.YLabel.String = '$dJ_n$';
sform = '$dJ_n = J_n - \frac{3(-1)^{1+n/2}}{(n+1)(n+3)}\left(\frac{l^2}{1+l^2}\right)^{n/2}$';
ah.Title.String = sform;
J_err = max(abs(dJ)) %#ok<NOPTS>