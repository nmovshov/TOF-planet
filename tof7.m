function [Js, out] = tof7(zvec, dvec, mrot, varargin)
%TOF7 Seventh-order Theory of Figures gravity coefficients.
%   Js = TOF7(zvec, dvec, mrot) returns 1-by-8 vector Js of gravity
%   coefficients J0 through J14 of a rotating fluid planet in hydrostatic
%   equilibrium. The mandatory inputs are a vector of mean radii zvec, vector
%   of corresponding densities dvec, and the rotation parameter mrot, assumed
%   normalized to the outer level surface mean radius.
%
%   [Js, out] = TOF7(zvec, dvec, mrot, 'NAME1',VALUE1, 'NAME2',VALUE2,...)
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
% Js : 1-by-8 vector, real
%     Even harmonic gravity coefficients J0 to J14 (J0 is included as a sanity
%     check and test of convergence).
% out : struct
%     A structure holding other quantities calculated in the course of running
%     tof, including the shape functions that define the full solution.
%
% Algorithm
% ---------
% Theory of figures equations and coefficients from Nettelmann et al., 2021.

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

% Load coefficients for ToF7 equations
try
    C7 = load(opts.C7file);
    C7 = C7.C7;
catch
    error('Failed to load coefficients from ''C7file''=%s',opts.C7file)
end

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
    ss.s10(N,1)=0; ss.s12(N,1)=0; ss.s14(N,1)=0;
else
    try
        ss = opts.ss_guesses;
        fs = B1617(ss, C7);
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
Js = [0, 0, 0, 0, 0, 0, 0, 0]; % J0=0 ensures at least one iteration
for iter=1:opts.maxiter
    % Equations B.16-B.17
    fs = skipnspline_B1617(ss, zvec, xind, C7);

    % Equation B.9
    SS = B9(zvec, dvec, fs);

    % And finally, the system of simultaneous equations B.12-B.15.
    ss = skipnspline_B1215(ss, SS, mrot, zvec, xind, C7);

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
    warning('TOF7:maxiter','Figure functions may not be fully converged.')
end

%% Return
Js = new_Js; % may as well use the latest...
out.dJs = dJs;
out.iter = iter;
out.a0 = a0;
out.qrot = mrot*a0^3;
out.ss = ss;
out.SS = SS;
out.A0 = B4(ss,SS,mrot,C7);

end

%% Helper functions
function print_usage()
fprintf('Usage:\n\ttof7(zvec,dvec,mrot,''name'',value)\n')
fprintf('Name-Value pair arguments:\n')
fprintf('tol - Convergence tolerance [ positive real {1e-10} ]\n')
fprintf('maxiter - Number of iterations allowed for relaxation to equilibrium shape [ positive integer {100} ]\n')
fprintf('xlevels - Solve shape functions on xlevels and spline the rest [ integer scalar or vector {-1} ]\n')
fprintf('ss_guesses - Initial guess for shape functions [ scalar struct {[]} ]\n')
end

function options = parsem(varargin)
p = inputParser;
p.FunctionName = 'tof7.m';

p.addParameter('C7file',mfilename('fullpath'))
p.addParameter('tol',1e-10,@(x)isscalar(x)&&isreal(x)&&x>0)
p.addParameter('maxiter',100,@(x)isscalar(x)&&isreal(x)&&x>0&&mod(x,1)==0)
p.addParameter('xlevels',-1,@(x)validateattributes(x,{'numeric'},{'vector','integer'}))
p.addParameter('ss_guesses',struct(),@(x)isscalar(x)&&isstruct(x))

% Parse name-value pairs and return
p.parse(varargin{:})
options = p.Results;

% Some post-processing required
if any(strcmp('C7file',p.UsingDefaults))
    [c7p,~,~] = fileparts(options.C7file);
    options.C7file = fullfile(c7p,'tof7_coeffs.mat');
end

end

function fs = B1617(ss, C7)
% The ToF7 equivalent of Nettelmann 2017 eqs. B.16 and B.17.

N = size(ss.s0,1);
sarray = [ss.s2, ss.s4, ss.s6, ss.s8, ss.s10, ss.s12, ss.s14];

fs.f0 = zeros(N,1);
block = C7.f.f0;
for k=1:size(block,1)
    if any(block(k,1:end-1))
        pows = repmat(block(k,1:end-1), [N,1]);
        c = block(k,end);
        fs.f0 = fs.f0 + c*prod((sarray.^pows),2);
    else
        fs.f0 = fs.f0 + block(k,end)*ones(N,1);
    end
end

fs.f2 = zeros(N,1);
block = C7.f.f2;
for k=1:size(block,1)
    if any(block(k,1:end-1))
        pows = repmat(block(k,1:end-1), [N,1]);
        c = block(k,end);
        fs.f2 = fs.f2 + c*prod((sarray.^pows),2);
    else
        fs.f2 = fs.f2 + block(k,end)*ones(N,1);
    end
end

fs.f4 = zeros(N,1);
block = C7.f.f4;
for k=1:size(block,1)
    if any(block(k,1:end-1))
        pows = repmat(block(k,1:end-1), [N,1]);
        c = block(k,end);
        fs.f4 = fs.f4 + c*prod((sarray.^pows),2);
    else
        fs.f4 = fs.f4 + block(k,end)*ones(N,1);
    end
end

fs.f6 = zeros(N,1);
block = C7.f.f6;
for k=1:size(block,1)
    if any(block(k,1:end-1))
        pows = repmat(block(k,1:end-1), [N,1]);
        c = block(k,end);
        fs.f6 = fs.f6 + c*prod((sarray.^pows),2);
    else
        fs.f6 = fs.f6 + block(k,end)*ones(N,1);
    end
end

fs.f8 = zeros(N,1);
block = C7.f.f8;
for k=1:size(block,1)
    if any(block(k,1:end-1))
        pows = repmat(block(k,1:end-1), [N,1]);
        c = block(k,end);
        fs.f8 = fs.f8 + c*prod((sarray.^pows),2);
    else
        fs.f8 = fs.f8 + block(k,end)*ones(N,1);
    end
end

fs.f10 = zeros(N,1);
block = C7.f.f10;
for k=1:size(block,1)
    if any(block(k,1:end-1))
        pows = repmat(block(k,1:end-1), [N,1]);
        c = block(k,end);
        fs.f10 = fs.f10 + c*prod((sarray.^pows),2);
    else
        fs.f10 = fs.f10 + block(k,end)*ones(N,1);
    end
end

fs.f12 = zeros(N,1);
block = C7.f.f12;
for k=1:size(block,1)
    if any(block(k,1:end-1))
        pows = repmat(block(k,1:end-1), [N,1]);
        c = block(k,end);
        fs.f12 = fs.f12 + c*prod((sarray.^pows),2);
    else
        fs.f12 = fs.f12 + block(k,end)*ones(N,1);
    end
end

fs.f14 = zeros(N,1);
block = C7.f.f14;
for k=1:size(block,1)
    if any(block(k,1:end-1))
        pows = repmat(block(k,1:end-1), [N,1]);
        c = block(k,end);
        fs.f14 = fs.f14 + c*prod((sarray.^pows),2);
    else
        fs.f14 = fs.f14 + block(k,end)*ones(N,1);
    end
end

fs.f0p = zeros(N,1);
block = C7.f.f0p;
for k=1:size(block,1)
    if any(block(k,1:end-1))
        pows = repmat(block(k,1:end-1), [N,1]);
        c = block(k,end);
        fs.f0p = fs.f0p + c*prod((sarray.^pows),2);
    else
        fs.f0p = fs.f0p + block(k,end)*ones(N,1);
    end
end

fs.f2p = zeros(N,1);
block = C7.f.f2p;
for k=1:size(block,1)
    if any(block(k,1:end-1))
        pows = repmat(block(k,1:end-1), [N,1]);
        c = block(k,end);
        fs.f2p = fs.f2p + c*prod((sarray.^pows),2);
    else
        fs.f2p = fs.f2p + block(k,end)*ones(N,1);
    end
end

fs.f4p = zeros(N,1);
block = C7.f.f4p;
for k=1:size(block,1)
    if any(block(k,1:end-1))
        pows = repmat(block(k,1:end-1), [N,1]);
        c = block(k,end);
        fs.f4p = fs.f4p + c*prod((sarray.^pows),2);
    else
        fs.f4p = fs.f4p + block(k,end)*ones(N,1);
    end
end

fs.f6p = zeros(N,1);
block = C7.f.f6p;
for k=1:size(block,1)
    if any(block(k,1:end-1))
        pows = repmat(block(k,1:end-1), [N,1]);
        c = block(k,end);
        fs.f6p = fs.f6p + c*prod((sarray.^pows),2);
    else
        fs.f6p = fs.f6p + block(k,end)*ones(N,1);
    end
end

fs.f8p = zeros(N,1);
block = C7.f.f8p;
for k=1:size(block,1)
    if any(block(k,1:end-1))
        pows = repmat(block(k,1:end-1), [N,1]);
        c = block(k,end);
        fs.f8p = fs.f8p + c*prod((sarray.^pows),2);
    else
        fs.f8p = fs.f8p + block(k,end)*ones(N,1);
    end
end

fs.f10p = zeros(N,1);
block = C7.f.f10p;
for k=1:size(block,1)
    if any(block(k,1:end-1))
        pows = repmat(block(k,1:end-1), [N,1]);
        c = block(k,end);
        fs.f10p = fs.f10p + c*prod((sarray.^pows),2);
    else
        fs.f10p = fs.f10p + block(k,end)*ones(N,1);
    end
end

fs.f12p = zeros(N,1);
block = C7.f.f12p;
for k=1:size(block,1)
    if any(block(k,1:end-1))
        pows = repmat(block(k,1:end-1), [N,1]);
        c = block(k,end);
        fs.f12p = fs.f12p + c*prod((sarray.^pows),2);
    else
        fs.f12p = fs.f12p + block(k,end)*ones(N,1);
    end
end

fs.f14p = zeros(N,1);
block = C7.f.f14p;
for k=1:size(block,1)
    if any(block(k,1:end-1))
        pows = repmat(block(k,1:end-1), [N,1]);
        c = block(k,end);
        fs.f14p = fs.f14p + c*prod((sarray.^pows),2);
    else
        fs.f14p = fs.f14p + block(k,end)*ones(N,1);
    end
end

end

function SS = B9(Z, D, fs)
% The ToF7 equivalent of Nettelmann 2017 eq. B.9.

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

I10 = cumtrapz(D, Z.^(10+3).*fs.f10);
SS.S10 = D.*fs.f10 - Z.^-(10+3).*I10;

I12 = cumtrapz(D, Z.^(12+3).*fs.f12);
SS.S12 = D.*fs.f12 - Z.^-(12+3).*I12;

I14 = cumtrapz(D, Z.^(14+3).*fs.f14);
SS.S14 = D.*fs.f14 - Z.^-(14+3).*I14;

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

I10p = cumtrapz(D, Z.^(2-10).*fs.f10p);
I10p = I10p(N) - I10p;
SS.S10p = -D.*fs.f10p + Z.^-(2-10).*(D(N)*fs.f10p(N) - I10p);

I12p = cumtrapz(D, Z.^(2-12).*fs.f12p);
I12p = I12p(N) - I12p;
SS.S12p = -D.*fs.f12p + Z.^-(2-12).*(D(N)*fs.f12p(N) - I12p);

I14p = cumtrapz(D, Z.^(2-14).*fs.f14p);
I14p = I14p(N) - I14p;
SS.S14p = -D.*fs.f14p + Z.^-(2-14).*(D(N)*fs.f14p(N) - I14p);

end

function fs = skipnspline_B1617(ss, zvec, xind, C7)
% Update the ToF7 equivalent of B.16-B.17 for new f_2n.

% Skip
zs.s0 = ss.s0(xind); zs.s2 = ss.s2(xind); zs.s4 = ss.s4(xind); zs.s6 = ss.s6(xind);
zs.s8 = ss.s8(xind); zs.s10 = ss.s10(xind); zs.s12 = ss.s12(xind); zs.s14 = ss.s14(xind);
fs = B1617(zs, C7);

% And spline
if length(xind) < length(zvec)
    fs.f0   = spline(zvec(xind), fs.f0, zvec);
    fs.f2   = spline(zvec(xind), fs.f2, zvec);
    fs.f4   = spline(zvec(xind), fs.f4, zvec);
    fs.f6   = spline(zvec(xind), fs.f6, zvec);
    fs.f8   = spline(zvec(xind), fs.f8, zvec);
    fs.f10  = spline(zvec(xind), fs.f10, zvec);
    fs.f12  = spline(zvec(xind), fs.f12, zvec);
    fs.f14  = spline(zvec(xind), fs.f14, zvec);
    fs.f0p  = spline(zvec(xind), fs.f0p, zvec);
    fs.f2p  = spline(zvec(xind), fs.f2p, zvec);
    fs.f4p  = spline(zvec(xind), fs.f4p, zvec);
    fs.f6p  = spline(zvec(xind), fs.f6p, zvec);
    fs.f8p  = spline(zvec(xind), fs.f8p, zvec);
    fs.f10p = spline(zvec(xind), fs.f10p, zvec);
    fs.f12p = spline(zvec(xind), fs.f12p, zvec);
    fs.f14p = spline(zvec(xind), fs.f14p, zvec);
end
end

function ss = skipnspline_B1215(ss0, SS, mrot, zvec, xind, C7)
% Update the ToF7 equivalent of the system B.12-B.15 for new s_2n.

% Skip
Zs = [SS.S0(xind), SS.S2(xind), SS.S4(xind), SS.S6(xind), SS.S8(xind),...
      SS.S10(xind), SS.S12(xind), SS.S14(xind)];
Zps = [SS.S0p(xind), SS.S2p(xind), SS.S4p(xind), SS.S6p(xind), SS.S8p(xind),...
       SS.S10p(xind), SS.S12p(xind), SS.S14p(xind)];
zs = [ss0.s2(xind), ss0.s4(xind), ss0.s6(xind), ss0.s8(xind),...
      ss0.s10(xind), ss0.s12(xind), ss0.s14(xind)];

newzs = B1215(zs, Zs, Zps, mrot, C7);
newz0 = A24(newzs);
Y = [newz0, newzs];

% And spline
if length(xind) < length(zvec)
    ss.s0  = spline(zvec(xind), Y(:,1), zvec);
    ss.s2  = spline(zvec(xind), Y(:,2), zvec);
    ss.s4  = spline(zvec(xind), Y(:,3), zvec);
    ss.s6  = spline(zvec(xind), Y(:,4), zvec);
    ss.s8  = spline(zvec(xind), Y(:,5), zvec);
    ss.s10 = spline(zvec(xind), Y(:,6), zvec);
    ss.s12 = spline(zvec(xind), Y(:,7), zvec);
    ss.s14 = spline(zvec(xind), Y(:,8), zvec);
else
    ss.s0  = Y(:,1);
    ss.s2  = Y(:,2);
    ss.s4  = Y(:,3);
    ss.s6  = Y(:,4);
    ss.s8  = Y(:,5);
    ss.s10 = Y(:,6);
    ss.s12 = Y(:,7);
    ss.s14 = Y(:,8);
end
end

function s0 = A24(sn)
% Nettelmann et al 2021 eq. A.24

s2 = sn(:,1); s4 = sn(:,2); s6 = sn(:,3);

s0 = 0*s2;
s0 = s0 - (1/5)*s2.^2;
s0 = s0 - (2/105)*s2.^3;
s0 = s0 - (1/9)*s4.^2 - (2/35)*s2.^2.*s4;
s0 = s0 - (2/525)*s2.^5 - (20/693)*s2.*s4.^2;
s0 = s0 + (127/55125)*s2.^6 - (2/175)*s2.^4.*s4 - (6/1001)*s4.^3 - ...
          (10/143)*s2.*s4.*s6 - (1/13)*s6.^2;
s0 = s0 - (2/2625)*s2.^7 - (82/10395)*s2.^3.*s4.^2 - (8/525)*s2.^5.*s4 - ...
          (14/715)*s2.*s6.^2 - (20/1257)*s4.^2.*s6;
end

function A0 = B4(ss, SS, m, C7)
S = [SS.S0, SS.S2, SS.S4, SS.S6, SS.S8, SS.S10, SS.S12, SS.S14];
Sp = [SS.S0p, SS.S2p, SS.S4p, SS.S6p, SS.S8p, SS.S10p, SS.S12p, SS.S14p];
s = [ss.s2, ss.s4, ss.s6, ss.s8, ss.s10, ss.s12, ss.s14];

mterm = (m/3)*ones(size(S,1),1);
terms = [S, Sp, mterm];
fields = fieldnames(C7.A.A0);
assert(size(terms,2) == length(fields))
nt = size(terms,2);
N = size(s,1);

A0 = 0*(1:N)';
for j=1:nt
    block = C7.A.A0.(fields{j});
    for k=1:size(block,1)
        if any(block(k,1:end-1))
            pows = repmat(block(k,1:end-1), N, 1);
            c = block(k,end);
            A0 = A0 + terms(:,j).*prod((s.^pows),2)*c;
        else
            A0 = A0 + terms(:,j)*block(k,end);
        end
    end
end
end

function new_s = B1215(s, S, Sp, m, C7)
% Compute ToF7 equivalent of the RHS of B.12-B.15 and "solve" for sn.

mterm = (m/3)*ones(size(S,1),1);
terms = [S, Sp, mterm];
fields = fieldnames(C7.A.A2); % same for all As
assert(size(terms,2) == length(fields))
nt = size(terms,2);
N = size(s,1);

% A2 (not including -s2S0)
A2 = 0*(1:N)';
for j=1:nt
    block = C7.A.A2.(fields{j});
    for k=1:size(block,1)
        if any(block(k,1:end-1))
            pows = repmat(block(k,1:end-1), N, 1);
            c = block(k,end);
            A2 = A2 + terms(:,j).*prod((s.^pows),2)*c;
        else
            A2 = A2 + terms(:,j)*block(k,end);
        end
    end
end

% A4 (not including -s4S0)
A4 = 0*(1:N)';
for j=1:nt
    block = C7.A.A4.(fields{j});
    for k=1:size(block,1)
        if any(block(k,1:end-1))
            pows = repmat(block(k,1:end-1), N, 1);
            c = block(k,end);
            A4 = A4 + terms(:,j).*prod((s.^pows),2)*c;
        else
            A4 = A4 + terms(:,j)*block(k,end);
        end
    end
end

% A6 (not including -s6S0)
A6 = 0*(1:N)';
for j=1:nt
    block = C7.A.A6.(fields{j});
    for k=1:size(block,1)
        if any(block(k,1:end-1))
            pows = repmat(block(k,1:end-1), N, 1);
            c = block(k,end);
            A6 = A6 + terms(:,j).*prod((s.^pows),2)*c;
        else
            A6 = A6 + terms(:,j)*block(k,end);
        end
    end
end

% A8 (not including -s8S0)
A8 = 0*(1:N)';
for j=1:nt
    block = C7.A.A8.(fields{j});
    for k=1:size(block,1)
        if any(block(k,1:end-1))
            pows = repmat(block(k,1:end-1), N, 1);
            c = block(k,end);
            A8 = A8 + terms(:,j).*prod((s.^pows),2)*c;
        else
            A8 = A8 + terms(:,j)*block(k,end);
        end
    end
end

% A10 (not including -s10S0)
A10 = 0*(1:N)';
for j=1:nt
    block = C7.A.A10.(fields{j});
    for k=1:size(block,1)
        if any(block(k,1:end-1))
            pows = repmat(block(k,1:end-1), N, 1);
            c = block(k,end);
            A10 = A10 + terms(:,j).*prod((s.^pows),2)*c;
        else
            A10 = A10 + terms(:,j)*block(k,end);
        end
    end
end

% A12 (not including -s12S0)
A12 = 0*(1:N)';
for j=1:nt
    block = C7.A.A12.(fields{j});
    for k=1:size(block,1)
        if any(block(k,1:end-1))
            pows = repmat(block(k,1:end-1), N, 1);
            c = block(k,end);
            A12 = A12 + terms(:,j).*prod((s.^pows),2)*c;
        else
            A12 = A12 + terms(:,j)*block(k,end);
        end
    end
end

% A14 (not including -s14S0)
A14 = 0*(1:N)';
for j=1:nt
    block = C7.A.A14.(fields{j});
    for k=1:size(block,1)
        if any(block(k,1:end-1))
            pows = repmat(block(k,1:end-1), N, 1);
            c = block(k,end);
            A14 = A14 + terms(:,j).*prod((s.^pows),2)*c;
        else
            A14 = A14 + terms(:,j)*block(k,end);
        end
    end
end

S0 = S(:,1);
new_s = [A2./S0, A4./S0, A6./S0, A8./S0, A10./S0, A12./S0, A14./S0];
end

function [Js, aos] = B111(ss, SS)
% Return Js from SS, with Req/Rm a necessary side effect.

N = length(ss.s0);
s0 = ss.s0(N); s2 = ss.s2(N); s4 = ss.s4(N); s6 = ss.s6(N); s8 = ss.s8(N);
s10 = ss.s10(N); s12 = ss.s12(N); s14 = ss.s14(N);
aos = 1 + s0 - (1/2)*s2 + (3/8)*s4 - (5/16)*s6 + (35/128)*s8 - ...
 (63/256)*s10 + (231/1024)*s12 - (429/2048)*s14;

J0  = -(aos^-0)*SS.S0(N);
J2  = -(aos^-2)*SS.S2(N);
J4  = -(aos^-4)*SS.S4(N);
J6  = -(aos^-6)*SS.S6(N);
J8  = -(aos^-8)*SS.S8(N);
J10 = -(aos^-10)*SS.S10(N);
J12 = -(aos^-12)*SS.S12(N);
J14 = -(aos^-14)*SS.S14(N);

Js = [J0, J2, J4, J6, J8, J10, J12, J14];
end
