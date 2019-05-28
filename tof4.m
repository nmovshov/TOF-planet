function [Js, out] = tof4(zvec, dvec, mrot, tol, maxiter, ss_guesses, sskip)
%TOF4 Forth-order Theory of Figures gravity coefficients.
%   Js = TOF4(zvec, dvec, mrot)
%   Js = TOF4(zvec, dvec, mrot, tol, maxiter, ss_guesses, sskip)

%% Input parsing
if nargin == 0 && nargout == 0
    help('tof4.m')
    return
end
narginchk(3,7);
if nargin < 4 || isempty(tol), tol = 1e-6; end
if nargin < 5 || isempty(maxiter), maxiter = 100; end
if nargin < 6 || isempty(ss_guesses), ss_guesses = struct(); end
if nargin < 7 || isempty(sskip), sskip = 0; end
validateattributes(zvec,{'numeric'},{'finite','nonnegative','vector'},'','zvec',1)
validateattributes(dvec,{'numeric'},{'finite','nonnegative','vector'},'','dvec',2)
validateattributes(mrot,{'numeric'},{'finite','nonnegative','scalar'},'','mrot',3)
validateattributes(tol,{'numeric'},{'finite','positive','scalar'},'','tol',4)
validateattributes(maxiter,{'numeric'},{'positive','scalar','integer'},'','maxiter',5)
validateattributes(ss_guesses,{'struct'},{'scalar'},'','ss_guesses',6)
validateattributes(sskip,{'numeric'},{'nonnegative','scalar','integer'},'','sskip',7)
assert(length(zvec) == length(dvec),...
    'length(zvec)=%d~=%d=length(dvec)',length(zvec),length(dvec))
[zvec, I] = sort(zvec);
dvec = dvec(I);
zvec = zvec(:); % now it's a column for sure
dvec = dvec(:); % now it's a column for sure
if zvec(1) == 0, zvec(1) = eps; end

%% Normalize radii and density
dro = [dvec(end); diff(flipud(dvec))];
m = sum(dro.*flipud(zvec).^3);
robar = m/zvec(end)^3;
zvec = zvec/zvec(end);
dvec = dvec/robar;

%% Initialize local variables
if nargin < 6 || isempty(fieldnames(ss_guesses))
    N = length(zvec);
    ss.s0(N,1)=0; ss.s2(N,1)=0; ss.s4(N,1)=0; ss.s6(N,1)=0; ss.s8(N,1)=0;
else
    try
        ss = ss_guesses;
        fs = B1617(ss);
        SS = B9(zvec, dvec, fs);
    catch ME
        error('Shape functions guess failed because:\n%s',ME.message)
    end
end

%% The loop, following Nettelmann (2017) Appendix B
Js = [0, 0, 0, 0, 0]; % J0=0 ensures at least one iteration
for iter=1:maxiter
    % Equations B.16-B.17
    fs = B1617(ss);

    % Equation B.9
    SS = B9(zvec, dvec, fs);

    % And finally, the system of simultaneous equations B.12-B.15.
    if sskip == 0
        ss = solve_B1215(ss, SS, mrot);
    else
        ss = skipnspline_B1215(ss, SS, mrot, zvec, sskip);
    end

    % The Js, by eqs. B.1 and B.11
    [new_Js, a0] = B111(ss, SS);

    % Check for convergence to terminate
    dJs = abs(Js - new_Js)./abs(Js+eps);
    if all(dJs < tol)
        break
    elseif iter < maxiter
        Js = new_Js;
    end
end
if iter == maxiter
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

N = length(Z); % x(N) is faster than x(end)

% m = nan(N,1);
% fun = @(z)interp1(Z, D.*Z.^2, z, 'pchip');
% for k=1:N
%     m(k) = integral(fun,0,Z(k));
% end
% m = 4*pi*m;
%
% m(1) = 4*pi/3*D(1)*Z(1)^3;
% for k=2:N
%     m(k) = m(k-1) + 4*pi/3*D(k)*(Z(k)^3 - Z(k-1)^3);
% end
%
% SS.S0 = m./(m(N)*Z.^3);

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

function ss = solve_B1215(ss0, SS, mrot)
% Solve the system B.12-B.15 for unknowns s2,s4,s6,s8.

opts = optimset();
opts.TolX = 1e-10;
opts.TolFun = 1e-10;
opts.Display = 'off';
Y = nan(length(SS.S0), 5);
Zs = [SS.S0, SS.S2, SS.S4, SS.S6, SS.S8];       % temp variable for parfor slicing
Zps = [SS.S0p, SS.S2p, SS.S4p, SS.S6p, SS.S8p]; % temp variable for parfor slicing
zs = [ss0.s2, ss0.s4, ss0.s6, ss0.s8];          % temp variable for parfor slicing
parfor k=1:length(SS.S0)
    S = Zs(k,:);
    Sp = Zps(k,:);
    s0 = zs(k,:);
    fun = @(x)B1215(x, S, Sp, mrot);
    XX = [0, fsolve(fun, s0, opts)];
    XX(1) = -1/5*XX(2)^2 - 2/105*XX(2)^3 - 1/9*XX(3)^2 - 2/35*XX(2)^2*XX(3);
    Y(k,:) = XX;
end
ss.s0 = Y(:,1); ss.s2 = Y(:,2); ss.s4 = Y(:,3); ss.s6 = Y(:,4); ss.s8 = Y(:,5);
end

function ss = skipnspline_B1215(ss0, SS, mrot, zvec, sskip)
% Solve the system B.12-B.15 for unknowns s2,s4,s6,s8.

opts = optimset();
opts.TolX = 1e-10;
opts.TolFun = 1e-10;
opts.Display = 'off';
N = length(SS.S0);
ind = 1:sskip:N;
Y = nan(length(ind), 5);
 % (temp variables for parfor slicing)
Zs = [SS.S0(ind), SS.S2(ind), SS.S4(ind), SS.S6(ind), SS.S8(ind)];
Zps = [SS.S0p(ind), SS.S2p(ind), SS.S4p(ind), SS.S6p(ind), SS.S8p(ind)];
zs = [ss0.s2(ind), ss0.s4(ind), ss0.s6(ind), ss0.s8(ind)];
parfor k=1:length(ind)
    S = Zs(k,:);
    Sp = Zps(k,:);
    s0 = zs(k,:);
    fun = @(x)B1215(x, S, Sp, mrot);
    XX = [0, fsolve(fun, s0, opts)];
    XX(1) = -1/5*XX(2)^2 - 2/105*XX(2)^3 - 1/9*XX(3)^2 - 2/35*XX(2)^2*XX(3);
    Y(k,:) = XX;
end
ss.s0 = spline(zvec(ind), Y(:,1), zvec);
ss.s2 = spline(zvec(ind), Y(:,2), zvec);
ss.s4 = spline(zvec(ind), Y(:,3), zvec);
ss.s6 = spline(zvec(ind), Y(:,4), zvec);
ss.s8 = spline(zvec(ind), Y(:,5), zvec);
end

function A = B1215(s, S, Sp, m)
% Compute the RHS of B.12-B.15.

s2 = s(1); s4 = s(2); s6 = s(3); s8 = s(4);
S0 = S(1); S2 = S(2); S4 = S(3); S6 = S(4); S8 = S(5);
S0p = Sp(1); S2p = Sp(2); S4p = Sp(3); S6p = Sp(4); S8p = Sp(5); %#ok<NASGU>

% B.12
A2 = 0;
A2 = A2 + S0*(-1*s2 + 2/7*s2^2 + 4/7*s2*s4 - 29/35*s2^3 + 100/693*s4^2 + ...
               454/1155*s2^4 - 36/77*s2^2*s4);
A2 = A2 + S2*(1 - 6/7*s2 - 6/7*s4 + 111/35*s2^2 - 1242/385*s2^3 + 144/77*s2*s4);
A2 = A2 + S4*(-10/7*s2 - 500/693*s4 + 180/77*s2^2);
A2 = A2 + S2p*(1 + 4/7*s2 + 1/35*s2^2 + 4/7*s4 - 16/105*s2^3 + 24/77*s2*s4);
A2 = A2 + S4p*(8/7*s2 + 72/77*s2^2 + 400/693*s4);
A2 = A2 + m/3*(-1 + 10/7*s2 + 9/35*s2^2 - 4/7*s4 + 20/77*s2*s4 - 26/105*s2^3);

% B.13
A4 = 0;
A4 = A4 + S0*(-1*s4 + 18/35*s2^2 - 108/385*s2^3 + 40/77*s2*s4 + ...
              90/143*s2*s6 + 162/1001*s4^2 + 16902/25025*s2^4 - ...
              7369/5005*s2^2*s4);
A4 = A4 + S2*(-54/35*s2 - 60/77*s4 + 648/385*s2^2 - 135/143*s6 + ...
              21468/5005*s2*s4 - 122688/25025*s2^3);
A4 = A4 + S4*(1 - 100/77*s2 - 810/1001*s4 + 6368/1001*s2^2);
A4 = A4 + S6*(-315/143*s2);
A4 = A4 + S2p*(36/35*s2 + 108/385*s2^2 + 40/77*s4 + 3578/5005*s2*s4 - ...
               36/175*s2^3 + 90/143*s6);
A4 = A4 + S4p*(1 + 80/77*s2 + 1346/1001*s2^2 + 648/1001*s4);
A4 = A4 + S6p*(270/143*s2);
A4 = A4 + m/3*(-36/35*s2 + 114/77*s4 + 18/77*s2^2 - 978/5005*s2*s4 + ...
               36/175*s2^3 - 90/143*s6);

% B.14
A6 = 0;
A6 = A6 + S0*(-s6 + 10/11*s2*s4 - 18/77*s2^3 + 28/55*s2*s6 + 72/385*s2^4 + ...
              20/99*s4^2 - 54/77*s2^2*s4);
A6 = A6 + S2*(-15/11*s4 + 108/77*s2^2 - 42/55*s6 - 144/77*s2^3 + 216/77*s2*s4);
A6 = A6 + S4*(-25/11*s2 - 100/99*s4 + 270/77*s2^2);
A6 = A6 + S6*(1 - 98/55*s2);
A6 = A6 + S2p*(10/11*s4 + 18/77*s2^2 + 36/77*s2*s4 + 28/55*s6);
A6 = A6 + S4p*(20/11*s2 + 108/77*s2^2 + 80/99*s4);
A6 = A6 + S6p*(1 + 84/55*s2);
A6 = A6 + m/3*(-10/11*s4 - 18/77*s2^2 + 34/77*s2*s4 + 82/55*s6);

% B.15
A8 = 0;
A8 = A8 + S0*(-1*s8 + 56/65*s2*s6 + 72/715*s2^4 + 490/1287*s4^2 - 84/143*s2^2*s4);
A8 = A8 + S2*(-84/65*s6 - 144/143*s2^3 + 336/143*s2*s4);
A8 = A8 + S4*(-2450/1287*s4 + 420/143*s2^2);
A8 = A8 + S6*(-196/65*s2);
A8 = A8 + S8*(1);
A8 = A8 + S2p*(56/65*s6 + 56/143*s2*s4);
A8 = A8 + S4p*(1960/1287*s4 + 168/143*s2^2);
A8 = A8 + S6p*(168/65*s2);
A8 = A8 + S8p*(1);
A8 = A8 + m/3*(-56/65*s6 - 56/143*s2*s4);

A = [A2, A4, A6, A8];
end

function [Js, aos] = B111(ss, SS)
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
