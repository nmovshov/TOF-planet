function [Js, out] = tof4(zvec, dvec, mrot, varargin)
%TOF4 Forth-order Theory of Figures gravity coefficients.
%   Js = TOF4(zvec, dvec, mrot) returns 1-by-5 vector Js of gravity
%   coefficients J0 through J8 of a rotating fluid planet in hydrostatic
%   equilibrium. The mandatory inputs are a vector of mean radii zvec, vector
%   of corresponding densities dvec, and the rotation parameter mrot, assumed
%   normalized to the outer level surface mean radius.
%
%   [Js, out] = TOF4(zvec, dvec, mrot, 'NAME1',VALUE1, 'NAME2',VALUE2,...)
%   accepts additional parametrs as NAME/VALUE pairs , and also returns an
%   output struct holding diagnostic values and additional derived quantities,
%   including the shape functions defining the full hydrostatic equlibirium
%   solution.
%
% Inputs, required
% ----------------
% zvec : 1d array, positive real
%     Mean radii of level surfaces where density is defined, indexed from the
%     outside in, i.e., zvec(1)=s0 is the radius of the outermost level
%     surface. Units of zvec are unimportant as values will be normalized to
%     outer radius.
% dvec : 1d array, positive real
%     Density on level surfaces. Units are unimportant as values will be
%     normalized to the mean (bulk) density. The density should be
%     monotonically non-increasing with zvec, but this is not enforced.
% mrot : scalar, nonnegative
%     Dimensionless rotation parameter. Recall m = w^2s0^3/GM.
%
% Inputs, NAME/VALUE pairs
% tol : scalar, positive, (tol=1e-10)
%     Convergence tolerance for absolute change in J2 in successive iterations
%     (but keep in mind truncation error of ToF).
% maxiter : scalar, positive, integer, (maxiter=100)
%     Maximum number of algorithm iterations.
% xlevels : scalar or vector, nonnegative, integer (xlevels=-1)
%     Levels whose shape will be explicitly calculated. The shape functions
%     will be explicitly calculated for these levels, and spline-interpolated
%     in between. This can result in significant speedup with minimal loss of
%     precision, if the xlevels are chosen by trial and error to fit the
%     required precision and the spacing of density levels. A scalar value is
%     interpreted as a number of xlevels to be uniformaly distributed among the
%     density levels. For example, a smooth-density 1024-level model can
%     benefit from almost 16x-speedup by specifying xlevels=64 while retaining
%     a 10^-6 relative precision on J2. A vector value is interpreted as
%     indices of levels to be used as xlevels. Skip-n-spline is recommended for
%     very high resolution density profiles, N>~10^4. Disable skip-n-spline by
%     passing xlevels=-1 rather than xlevels=N, to avoid the spline ovearhead.
% ss_guesses : struct (default empty)
%     Initial guess for shape functions. This is not all that helpful in
%     speeding up convergence. It's occasionally helpful to preserve state
%     between successive calls.
%
% Outputs
% -------
% Js : 1-by-5 vector, real
%     Even harmonic gravity coefficients J0 to J8 (J0 is included as a sanity
%     check and test of convergence).
% out : struct
%     A structure holding other quantities calculated in the course of running
%     tof, including the shape functions that define the full solution.
%
% Algorithm
% ---------
% Theory of figures equations and coefficients from Nettelmann 2017 Appendix B.

%% Input parsing
% Zero inputs case, usage only
if nargin == 0
    print_usage()
    return
end
narginchk(3,inf);

% Mandatory inputs
validateattributes(zvec,{'numeric'},{'finite','nonnegative','vector'},'','zvec',1)
validateattributes(dvec,{'numeric'},{'finite','nonnegative','vector'},'','dvec',2)
validateattributes(mrot,{'numeric'},{'finite','nonnegative','scalar'},'','mrot',3)
assert(length(zvec) == length(dvec),...
    'length(zvec)=%d~=%d=length(dvec)',length(zvec),length(dvec))
[zvec, I] = sort(zvec);
dvec = dvec(I);
zvec = zvec(:); % now it's a column for sure
dvec = dvec(:); % now it's a column for sure
if zvec(1) == 0, zvec(1) = eps; end

% Optional arguments
opts = parsem(varargin{:});

% Normalize radii and density
dro = [dvec(end); diff(flipud(dvec))];
m = sum(dro.*flipud(zvec).^3);
robar = m/zvec(end)^3;
zvec = zvec/zvec(end);
dvec = dvec/robar;

%% Define and initialize local variables
nlay = length(zvec);

if isempty(fieldnames(opts.ss_guesses))
    N = nlay;
    ss.s0(N,1)=0; ss.s2(N,1)=0; ss.s4(N,1)=0; ss.s6(N,1)=0; ss.s8(N,1)=0;
else
    try
        ss = opts.ss_guesses;
        fs = B1617(ss);
        SS = B9(zvec, dvec, fs);
    catch ME
        error('Shape functions guess failed because:\n%s',ME.message)
    end
end

% Define down-sampled variabels (for skip-n-spline)
if isscalar(opts.xlevels)
    sskip = max(fix(nlay/opts.xlevels), 1);
    xind = 1:sskip:nlay;
else
    warning('Experimental feature, use with care.')
    xind = opts.xlevels;
end

%% The loop, following Nettelmann (2017) Appendix B
Js = [0, 0, 0, 0, 0]; % J0=0 ensures at least one iteration
for iter=1:opts.maxiter
    % Equations B.16-B.17
    fs = B1617(ss);

    % Equation B.9
    SS = B9(zvec, dvec, fs);

    % And finally, the system of simultaneous equations B.12-B.15.
    ss = skipnspline_B1215(ss, SS, mrot, zvec, xind);

    % The Js, by eqs. B.1 and B.11
    [new_Js, a0] = B111(ss, SS);

    % Check for convergence to terminate
    dJs = abs(Js - new_Js);
    if (iter > 1) && (dJs(2) < opts.tol)
        break
    elseif iter < opts.maxiter
        Js = new_Js;
    end
end
if iter == opts.maxiter
    warning('TOF4:maxiter','Figure functions may not be fully converged.')
end

%% Return
Js = new_Js; % may as well use the latest...
out.dJs = dJs;
out.iter = iter;
out.a0 = a0;
out.qrot = mrot*a0^3;
out.ss = ss;
out.SS = SS;
out.A0 = B4(ss,SS,mrot);

end

%% Helper functions
function print_usage()
fprintf('Usage:\n\ttof4(zvec,dvec,mrot,''name'',value)\n')
fprintf('Name-Value pair arguments:\n')
fprintf('tol - Convergence tolerance [ positive real {1e-10} ]\n')
fprintf('maxiter - Number of iterations allowed for relaxation to equilibrium shape [ positive integer {100} ]\n')
fprintf('xlevels - Solve shape functions on xlevels and spline the rest [ integer scalar or vector {-1} ]\n')
fprintf('ss_guesses - Initial guess for shape functions [ scalar struct {[]} ]\n')
end

function options = parsem(varargin)
p = inputParser;
p.FunctionName = 'tof4.m';

p.addParameter('tol',1e-10,@(x)isscalar(x)&&isreal(x)&&x>0)
p.addParameter('maxiter',100,@(x)isscalar(x)&&isreal(x)&&x>0&&mod(x,1)==0)
p.addParameter('xlevels',-1,@(x)validateattributes(x,{'numeric'},{'vector','integer'}))
p.addParameter('ss_guesses',struct(),@(x)isscalar(x)&&isstruct(x))

% Parse name-value pairs and return
p.parse(varargin{:})
options = p.Results;
end

function fs = B1617(ss)
% Nettelmann 2017 eqs. B.16 and B.17.

s0 = ss.s0; s2 = ss.s2; s4 = ss.s4; s6 = ss.s6; s8 = ss.s8; %#ok<NASGU>

fs.f0 = ones(size(ss.s0));

fs.f2 = (3/5)*s2 + (12/35)*s2.^2 + (6/175)*s2.^3 + (24/35)*s2.*s4 + ...
        (40/231)*s4.^2 + (216/385)*s2.^2.*s4 - (184/1925)*s2.^4;

fs.f4 = (1/3)*s4 + (18/35)*s2.^2 + (40/77)*s2.*s4 + (36/77)*s2.^3 + ...
        (90/143)*s2.*s6 + (162/1001)*s4.^2 + (6943/5005)*s2.^2.*s4 + ...
        (486/5005)*s2.^4;

fs.f6 = (3/13)*s6 + (120/143)*s2.*s4 + (72/143)*s2.^3 + (336/715)*s2.*s6 + ...
        (80/429)*s4.^2 + (216/143)*s2.^2.*s4 + (432/715)*s2.^4;

fs.f8 = (3/17)*s8 + (168/221)*s2.*s6 + (2450/7293)*s4.^2 + ...
        (3780/2431)*s2.^2.*s4 + (1296/2431)*s2.^4;

fs.f0p = (3/2) - (3/10)*s2.^2 - (2/35)*s2.^3 - (1/6)*s4.^2 - ...
         (6/35)*s2.^2.*s4 + (3/50)*s2.^4;

fs.f2p = (3/5)*s2 - (3/35)*s2.^2 - (6/35)*s2.*s4 + (36/175)*s2.^3 - ...
         (10/231)*s4.^2 - (17/275)*s2.^4 + (36/385)*s2.^2.*s4;

fs.f4p = (1/3)*s4 - (9/35)*s2.^2 - (20/77)*s2.*s4 - (45/143)*s2.*s6 - ...
         (81/1001)*s4.^2 + (1/5)*s2.^2.*s4;

fs.f6p = (3/13)*s6 - (75/143)*s2.*s4 + (270/1001)*s2.^3 - (50/429)*s4.^2 + ...
         (810/1001)*s2.^2.*s4 - (54/143)*s2.^4 - (42/143)*s2.*s6;

fs.f8p = (3/17)*s8 - (588/1105)*s2.*s6 - (1715/7293)*s4.^2 + ...
         (2352/2431)*s2.^2.*s4 - (4536/12155)*s2.^4;
end

function SS = B9(Z, D, fs)
% Nettelmann 2017 eq. B.9.

% TODO: replace cumtrapz call with cumtrapz one-liner

N = length(Z); % x(N) is faster than x(end)

I0 = cumtrapz(D, Z.^(0+3).*fs.f0);
SS.S0 = D.*fs.f0 - Z.^-(0+3).*I0;

I2 = cumtrapz(D, Z.^(2+3).*fs.f2);
SS.S2 = D.*fs.f2 - Z.^-(2+3).*I2;

I4 = cumtrapz(D, Z.^(4+3).*fs.f4);
SS.S4 = D.*fs.f4 - Z.^-(4+3).*I4;

I6 = cumtrapz(D, Z.^(6+3).*fs.f6);
SS.S6 = D.*fs.f6 - Z.^-(6+3).*I6;

I8 = cumtrapz(D, Z.^(8+3).*fs.f8);
SS.S8 = D.*fs.f8 - Z.^-(8+3).*I8;

I0p = cumtrapz(D, Z.^(2-0).*fs.f0p);
I0p = I0p(N) - I0p;
SS.S0p = -D.*fs.f0p + Z.^-(2-0).*(D(N)*fs.f0p(N) - I0p);

I2p = cumtrapz(D, Z.^(2-2).*fs.f2p);
I2p = I2p(N) - I2p;
SS.S2p = -D.*fs.f2p + Z.^-(2-2).*(D(N)*fs.f2p(N) - I2p);

I4p = cumtrapz(D, Z.^(2-4).*fs.f4p);
I4p = I4p(N) - I4p;
SS.S4p = -D.*fs.f4p + Z.^-(2-4).*(D(N)*fs.f4p(N) - I4p);

I6p = cumtrapz(D, Z.^(2-6).*fs.f6p);
I6p = I6p(N) - I6p;
SS.S6p = -D.*fs.f6p + Z.^-(2-6).*(D(N)*fs.f6p(N) - I6p);

I8p = cumtrapz(D, Z.^(2-8).*fs.f8p);
I8p = I8p(N) - I8p;
SS.S8p = -D.*fs.f8p + Z.^-(2-8).*(D(N)*fs.f8p(N) - I8p);

end

function A0 = B4(ss, SS, mrot)
% Compute the RHS of B.4 in Nettelmann (2017).

s2 = ss.s2; s4 = ss.s4;
S0 = SS.S0; S2 = SS.S2; S4 = SS.S4;
S0p = SS.S0p; S2p = SS.S2p; S4p = SS.S4p;

A0 = zeros(size(s2));
A0 = A0 + S0.*(1 + (2/5)*s2.^2 - (4/105)*s2.^3 + (2/9)*s4.^2 + (43/175)*s2.^4 - (4/35)*s2.^2.*s4);
A0 = A0 + S2.*(-(3/5)*s2 + (12/35)*s2.^2 - (234/175)*s2.^3 + (24/35)*s2.*s4);
A0 = A0 + S4.*((6/7)*s2.^2 - (5/9)*s4);
A0 = A0 + S0p.*(1);
A0 = A0 + S2p.*((2/5)*s2 + (2/35)*s2.^2 + (4/35)*s2.*s4 - (2/25)*s2.^3);
A0 = A0 + S4p.*((4/9)*s4 + (12/35)*s2.^2);
A0 = A0 + (mrot/3)*(1 - (2/5)*s2 - (9/35)*s2.^2 - (4/35)*s2.*s4 + (22/525)*s2.^3);
end

function ss = skipnspline_B1215(ss0, SS, mrot, zvec, xind)
% Update the system B.12-B.15 for new s2,s4,s6,s8.

% Skip
Zs = [SS.S0(xind), SS.S2(xind), SS.S4(xind), SS.S6(xind), SS.S8(xind)];
Zps = [SS.S0p(xind), SS.S2p(xind), SS.S4p(xind), SS.S6p(xind), SS.S8p(xind)];
zs = [ss0.s2(xind), ss0.s4(xind), ss0.s6(xind), ss0.s8(xind)];

newzs = B1215(zs, Zs, Zps, mrot);
newz0 = (-1/5)*newzs(:,1).^2 - (2/105)*newzs(:,1).^3 - ...
        (1/9)*newzs(:,2).^2 - 2/35*newzs(:,1).^2.*newzs(:,2);
Y = [newz0, newzs];

% And spline
if length(xind) < length(zvec)
    ss.s0 = spline(zvec(xind), Y(:,1), zvec);
    ss.s2 = spline(zvec(xind), Y(:,2), zvec);
    ss.s4 = spline(zvec(xind), Y(:,3), zvec);
    ss.s6 = spline(zvec(xind), Y(:,4), zvec);
    ss.s8 = spline(zvec(xind), Y(:,5), zvec);
else
    ss.s0 = Y(:,1);
    ss.s2 = Y(:,2);
    ss.s4 = Y(:,3);
    ss.s6 = Y(:,4);
    ss.s8 = Y(:,5);
end
end

function new_s = B1215(s, S, Sp, m)
% Compute the RHS of B.12-B.15 and "solve" for sn.

s2 = s(:,1); s4 = s(:,2); s6 = s(:,3); s8 = s(:,4);
S0 = S(:,1); S2 = S(:,2); S4 = S(:,3); S6 = S(:,4); S8 = S(:,5);
S0p = Sp(:,1); S2p = Sp(:,2); S4p = Sp(:,3); S6p = Sp(:,4); S8p = Sp(:,5);

% B.12 (not including -s2S0)
A2 = 0;
A2 = A2 + S0.*(2/7*s2.^2 + 4/7*s2.*s4 - 29/35*s2.^3 + 100/693*s4.^2 + ...
               454/1155*s2.^4 - 36/77*s2.^2.*s4);
A2 = A2 + S2.*(1 - 6/7*s2 - 6/7*s4 + 111/35*s2.^2 - 1242/385*s2.^3 + 144/77*s2.*s4);
A2 = A2 + S4.*(-10/7*s2 - 500/693*s4 + 180/77*s2.^2);
A2 = A2 + S2p.*(1 + 4/7*s2 + 1/35*s2.^2 + 4/7*s4 - 16/105*s2.^3 + 24/77*s2.*s4);
A2 = A2 + S4p.*(8/7*s2 + 72/77*s2.^2 + 400/693*s4);
A2 = A2 + m/3*(-1 + 10/7*s2 + 9/35*s2.^2 - 4/7*s4 + 20/77*s2.*s4 - 26/105*s2.^3);

% B.13 (not including -s4S0)
A4 = 0;
A4 = A4 + S0.*(18/35*s2.^2 - 108/385*s2.^3 + 40/77*s2.*s4 + ...
              90/143*s2.*s6 + 162/1001*s4.^2 + 16902/25025*s2.^4 - ...
              7369/5005*s2.^2.*s4);
A4 = A4 + S2.*(-54/35*s2 - 60/77*s4 + 648/385*s2.^2 - 135/143*s6 + ...
              21468/5005*s2.*s4 - 122688/25025*s2.^3);
A4 = A4 + S4.*(1 - 100/77*s2 - 810/1001*s4 + 6368/1001*s2.^2);
A4 = A4 + S6.*(-315/143*s2);
A4 = A4 + S2p.*(36/35*s2 + 108/385*s2.^2 + 40/77*s4 + 3578/5005*s2.*s4 - ...
               36/175*s2.^3 + 90/143*s6);
A4 = A4 + S4p.*(1 + 80/77*s2 + 1346/1001*s2.^2 + 648/1001*s4);
A4 = A4 + S6p.*(270/143*s2);
A4 = A4 + m/3*(-36/35*s2 + 114/77*s4 + 18/77*s2.^2 - 978/5005*s2.*s4 + ...
               36/175*s2.^3 - 90/143*s6);

% B.14 (not including -s6S0)
A6 = 0;
A6 = A6 + S0.*(10/11*s2.*s4 - 18/77*s2.^3 + 28/55*s2.*s6 + 72/385*s2.^4 + ...
              20/99*s4.^2 - 54/77*s2.^2.*s4);
A6 = A6 + S2.*(-15/11*s4 + 108/77*s2.^2 - 42/55*s6 - 144/77*s2.^3 + 216/77*s2.*s4);
A6 = A6 + S4.*(-25/11*s2 - 100/99*s4 + 270/77*s2.^2);
A6 = A6 + S6.*(1 - 98/55*s2);
A6 = A6 + S2p.*(10/11*s4 + 18/77*s2.^2 + 36/77*s2.*s4 + 28/55*s6);
A6 = A6 + S4p.*(20/11*s2 + 108/77*s2.^2 + 80/99*s4);
A6 = A6 + S6p.*(1 + 84/55*s2);
A6 = A6 + m/3*(-10/11*s4 - 18/77*s2.^2 + 34/77*s2.*s4 + 82/55*s6);

% B.15 (not including -s8S0)
A8 = 0;
A8 = A8 + S0.*(56/65*s2.*s6 + 72/715*s2.^4 + 490/1287*s4.^2 - 84/143*s2.^2.*s4);
A8 = A8 + S2.*(-84/65*s6 - 144/143*s2.^3 + 336/143*s2.*s4);
A8 = A8 + S4.*(-2450/1287*s4 + 420/143*s2.^2);
A8 = A8 + S6.*(-196/65*s2);
A8 = A8 + S8*(1);
A8 = A8 + S2p.*(56/65*s6 + 56/143*s2.*s4);
A8 = A8 + S4p.*(1960/1287*s4 + 168/143*s2.^2);
A8 = A8 + S6p.*(168/65*s2);
A8 = A8 + S8p*(1);
A8 = A8 + m/3*(-56/65*s6 - 56/143*s2.*s4);

new_s = [A2./S0, A4./S0, A6./S0, A8./S0];
end

function [Js, aos] = B111(ss, SS)
% Return Js from SS, with Req/Rm a necessary bonus.

N = length(ss.s0);
s0 = ss.s0(N); s2 = ss.s2(N); s4 = ss.s4(N); s6 = ss.s6(N); s8 = ss.s8(N);
aos = 1 + s0 - (1/2)*s2 + (3/8)*s4 - (5/16)*s6 + (35/128)*s8;

J0 = -(aos^-0)*SS.S0(N);
J2 = -(aos^-2)*SS.S2(N);
J4 = -(aos^-4)*SS.S4(N);
J6 = -(aos^-6)*SS.S6(N);
J8 = -(aos^-8)*SS.S8(N);

Js = [J0, J2, J4, J6, J8];
end
