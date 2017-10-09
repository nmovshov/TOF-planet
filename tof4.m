function [Js, out] = tof4(zvec, dvec, mrot, tol, maxiter)
%TOF4 Forth-order Theory of Figures gravity coefficients.
%   Js = TOF4(zvec, dvec, mrot)
%   Js = TOF4(zvec, dvec, mrot, tol, maxiter)

%% Input parsing
try
    narginchk(3,5);
catch ME
    help('tof4.m')
    rethrow(ME)
end
if nargin < 4 || isempty(tol), tol = 1e-6; end
if nargin < 5 || isempty(maxiter), maxiter = 100; end
validateattributes(zvec,{'numeric'},{'finite','nonnegative','vector'},'','zvec',1)
validateattributes(dvec,{'numeric'},{'finite','nonnegative','vector'},'','dvec',2)
validateattributes(mrot,{'numeric'},{'finite','nonnegative','scalar'},'','mrot',3)
validateattributes(tol,{'numeric'},{'finite','positive','scalar'},'','tol',4)
assert(length(zvec) == length(dvec),...
    'length(zvec)=%d~=%d=length(dvec)',length(zvec),length(dvec))
[zvec, I] = sort(zvec);
dvec = dvec(I);
assert(zvec(end) == 1,['Are you using dimensioned variables? ',...
    'Normalize radii to the outer radius and densities to the mean density.']);
zvec = zvec(:); % now it's a column for sure
dvec = dvec(:); % now it's a column for sure
if zvec(1) == 0, zvec(1) = eps; end

%% Initialize local variables
N = length(zvec);
ss.s0(N,1)=0; ss.s2(N,1)=0; ss.s4(N,1)=0; ss.s6(N,1)=0; ss.s8(N,1)=0;
Js = [0, 0, 0, 0, 0];

%% The loop, following Nettelmann (2017) Appendix B
for iter=1:maxiter
    % Equations B.16-B.17
    fs = B1617(ss);
    
    % Equation B.9
    SS = B9(zvec, dvec, fs);
    
    % And finally, the system of simultaneous equations B.12-B.15.
    ss = solve_B1215(ss, SS, mrot);
    
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
    warning('Figure functions may not be fully converged.')
end

%% Return
Js = new_Js; % may as well use the latest...
out.dJs = dJs;
out.iter = iter;
out.a0 = a0;
out.qrot = mrot*a0^3;
out.ss = ss;

end

function fs = B1617(ss)
% Nettelmann 2017 eqs. B.16 and B.17.
s0 = ss.s0; s2 = ss.s2; s4 = ss.s4; s6 = ss.s6; s8 = ss.s8;

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
         (81/1001)*s4.*2 + (1/5)*s2.^2.*s4;

fs.f6p = (3/13)*s6 - (75/143)*s2.*s4 + (270/1001)*s2.^3 - (50/429)*s4.^2 + ...
         (810/1001)*s2.^2.*s4 - (54/143)*s2.^4 - (42/143)*s2.*s6;

fs.f8p = (3/17)*s8 - (588/1105)*s2.*s6 - (1715/7293)*s4.^2 + ...
         (2352/2431)*s2.^2.*s4 - (4536/12155)*s2.^4;
end

function SS = B9(Z, D, fs)
% Nettelmann 2017 eq. B.9.

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

function ss = solve_B1215(ss0, SS, mrot)
% Solve the system B.12-B.15 for unknowns s2,s4,s6,s8.

opts = optimset();
opts.display = 'off';
Y = nan(length(SS.S0), 5);
for k=1:length(SS.S0)
    S = [SS.S0(k), SS.S2(k), SS.S4(k), SS.S6(k), SS.S8(k)];
    Sp = [SS.S0p(k), SS.S2p(k), SS.S4p(k), SS.S6p(k), SS.S8p(k)];
    s0 = [ss0.s2(k), ss0.s4(k), ss0.s6(k), ss0.s8(k)];
    fun = @(x)B1215(x, S, Sp, mrot);
    Y(k,2:5) = fsolve(fun, s0, opts);

    % And eq. B.2 too
    Y(k,1) = -1/5*Y(k,2)^2 - 2/105*Y(k,2)^3 - 1/9*Y(k,3)^2 - 2/35*Y(k,2)^2*Y(k,3);
end
ss.s0 = Y(:,1); ss.s2 = Y(:,2); ss.s4 = Y(:,3); ss.s6 = Y(:,4); ss.s8 = Y(:,5);
end

function A = B1215(s, S, Sp, m)
% Compute the RHS of B.12-B.15.

s2 = s(1); s4 = s(2); s6 = s(3); s8 = s(4);
S0 = S(1); S2 = S(2); S4 = S(3); S6 = S(4); S8 = S(5);
S0p = Sp(1); S2p = Sp(2); S4p = Sp(3); S6p = Sp(4); S8p = Sp(5);

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

%% External functions
% Below I include customized versions of built-in MATLAB functions. This ensures
% that these optimized versions will be used in calls to tof4 no matter what else
% is on the MATLAB path.
%
% * fsolve() - commenting out creation of detailed exit messages to reduce
%              overhead.
% * trustnleqn(), dogleg() - dependencies of fsolve(), unmodifed.
% * optimget() - optimized version of optimget with zero error checks!
% * optimset() - set all values for fsolve

function defaultopt = optimset()
defaultopt = struct(...
    'Algorithm','trust-region-dogleg',...
    'DerivativeCheck','off',...
    'Diagnostics','off',...
    'DiffMaxChange',Inf,...
    'DiffMinChange',0,...
    'Display','final',...
    'FinDiffRelStep', [], ...
    'FinDiffType','forward',...
    'FunValCheck','off',...
    'InitDamping', 0.01, ...
    'Jacobian','off',...
    'JacobMult',[],... 
    'JacobPattern','sparse(ones(Jrows,Jcols))',...
    'MaxFunEvals',[],...
    'MaxIter',400,...
    'MaxPCGIter','max(1,floor(numberOfVariables/2))',...
    'OutputFcn',[],...
    'PlotFcns',[],...
    'PrecondBandWidth',Inf,...
    'ScaleProblem','none',...
    'TolFun', 1e-6,...
    'TolFunValue', 1e-6, ...    
    'TolPCG',0.1,...
    'TolX',1e-6,...
    'TypicalX','ones(numberOfVariables,1)', ...
    'UseParallel', false );
end

function o = optimget(options,name,~,~)
o = options.(name);
end

function [x,FVAL,EXITFLAG,OUTPUT,JACOB] = fsolve(FUN,x,options,varargin)
%FSOLVE solves systems of nonlinear equations of several variables.
%
%   FSOLVE attempts to solve equations of the form:
%             
%   F(X) = 0    where F and X may be vectors or matrices.   
%
%   FSOLVE implements three different algorithms: trust region dogleg,
%   trust region, and Levenberg-Marquardt. Choose one via the option
%   Algorithm: for instance, to choose trust region, set OPTIONS =
%   optimoptions('fsolve','Algorithm','trust-region'), and then pass
%   OPTIONS to FSOLVE.
%    
%   X = FSOLVE(FUN,X0) starts at the matrix X0 and tries to solve the 
%   equations in FUN.  FUN accepts input X and returns a vector (matrix) of 
%   equation values F evaluated at X. 
%
%   X = FSOLVE(FUN,X0,OPTIONS) solves the equations with the default
%   optimization parameters replaced by values in OPTIONS, an argument
%   created with the OPTIMOPTIONS function.  See OPTIMOPTIONS for details.
%   Use the SpecifyObjectiveGradient option to specify that FUN also
%   returns a second output argument J that is the Jacobian matrix at the
%   point X. If FUN returns a vector F of m components when X has length n,
%   then J is an m-by-n matrix where J(i,j) is the partial derivative of
%   F(i) with respect to x(j). (Note that the Jacobian J is the transpose
%   of the gradient of F.)
%
%   X = FSOLVE(PROBLEM) solves system defined in PROBLEM. PROBLEM is a
%   structure with the function FUN in PROBLEM.objective, the start point
%   in PROBLEM.x0, the options structure in PROBLEM.options, and solver
%   name 'fsolve' in PROBLEM.solver.  Use this syntax to solve at the 
%   command line a problem exported from OPTIMTOOL. 
%
%   [X,FVAL] = FSOLVE(FUN,X0,...) returns the value of the equations FUN 
%   at X. 
%
%   [X,FVAL,EXITFLAG] = FSOLVE(FUN,X0,...) returns an EXITFLAG that
%   describes the exit condition. Possible values of EXITFLAG and the
%   corresponding exit conditions are listed below. See the documentation
%   for a complete description.
%
%     1  FSOLVE converged to a root.
%     2  Change in X too small.
%     3  Change in residual norm too small.
%     4  Computed search direction too small.
%     0  Too many function evaluations or iterations.
%    -1  Stopped by output/plot function.
%    -2  Converged to a point that is not a root.
%    -3  Trust region radius too small (Trust-region-dogleg).
%
%   [X,FVAL,EXITFLAG,OUTPUT] = FSOLVE(FUN,X0,...) returns a structure 
%   OUTPUT with the number of iterations taken in OUTPUT.iterations, the 
%   number of function evaluations in OUTPUT.funcCount, the algorithm used 
%   in OUTPUT.algorithm, the number of CG iterations (if used) in 
%   OUTPUT.cgiterations, the first-order optimality (if used) in 
%   OUTPUT.firstorderopt, and the exit message in OUTPUT.message.
%
%   [X,FVAL,EXITFLAG,OUTPUT,JACOB] = FSOLVE(FUN,X0,...) returns the 
%   Jacobian of FUN at X.  
%
%   Examples
%     FUN can be specified using @:
%        x = fsolve(@myfun,[2 3 4],optimoptions('fsolve','Display','iter'))
%
%   where myfun is a MATLAB function such as:
%
%       function F = myfun(x)
%       F = sin(x);
%
%   FUN can also be an anonymous function:
%
%       x = fsolve(@(x) sin(3*x),[1 4],optimoptions('fsolve','Display','off'))
%
%   If FUN is parameterized, you can use anonymous functions to capture the 
%   problem-dependent parameters. Suppose you want to solve the system of 
%   nonlinear equations given in the function myfun, which is parameterized 
%   by its second argument c. Here myfun is a MATLAB file function such as
%     
%       function F = myfun(x,c)
%       F = [ 2*x(1) - x(2) - exp(c*x(1))
%             -x(1) + 2*x(2) - exp(c*x(2))];
%           
%   To solve the system of equations for a specific value of c, first 
%   assign the value to c. Then create a one-argument anonymous function 
%   that captures that value of c and calls myfun with two arguments. 
%   Finally, pass this anonymous function to FSOLVE:
%
%       c = -1; % define parameter first
%       x = fsolve(@(x) myfun(x,c),[-5;-5])
%
%   See also OPTIMOPTIONS, LSQNONLIN, @.

%   Copyright 1990-2015 The MathWorks, Inc.

% ------------Initialization----------------
defaultopt = struct(...
    'Algorithm','trust-region-dogleg',...
    'DerivativeCheck','off',...
    'Diagnostics','off',...
    'DiffMaxChange',Inf,...
    'DiffMinChange',0,...
    'Display','final',...
    'FinDiffRelStep', [], ...
    'FinDiffType','forward',...
    'FunValCheck','off',...
    'InitDamping', 0.01, ...
    'Jacobian','off',...
    'JacobMult',[],... 
    'JacobPattern','sparse(ones(Jrows,Jcols))',...
    'MaxFunEvals',[],...
    'MaxIter',400,...
    'MaxPCGIter','max(1,floor(numberOfVariables/2))',...
    'OutputFcn',[],...
    'PlotFcns',[],...
    'PrecondBandWidth',Inf,...
    'ScaleProblem','none',...
    'TolFun', 1e-6,...
    'TolFunValue', 1e-6, ...    
    'TolPCG',0.1,...
    'TolX',1e-6,...
    'TypicalX','ones(numberOfVariables,1)', ...
    'UseParallel', false );

% If just 'defaults' passed in, return the default options in X
if nargin == 1 && nargout <= 1 && strcmpi(FUN,'defaults')
   x = defaultopt;
   return
end

if nargin < 3, options=[]; end

% Detect problem structure input
if nargin == 1
    if isa(FUN,'struct')
        [FUN,x,options] = separateOptimStruct(FUN);
    else % Single input and non-structure.
        error(message('optim:fsolve:InputArg'));
    end
end

% The 'trust-region-reflective' algorithm has been renamed to
% 'trust-region'. Options objects will handle this conversion. However, we
% need to check any structures which may undergo no checking before this
% point.
if isstruct(options) && isfield(options, 'Algorithm') && ...
        ~iscell(options.Algorithm) && ...
        strcmp(options.Algorithm, 'trust-region-reflective')
    options.Algorithm = 'trust-region';
end

% Prepare the options for the solver
[options, optionFeedback] = prepareOptionsForSolver(options, 'fsolve');

if nargin == 0
  error(message('optim:fsolve:NotEnoughInputs'))
end
% Check for non-double inputs
msg = isoptimargdbl('FSOLVE', {'X0'}, x);
if ~isempty(msg)
    error('optim:fsolve:NonDoubleInput',msg);
end

LB = []; UB = [];
sizes.xShape = size(x);
xstart = x(:);
sizes.nVar = length(xstart);
sizes.mNonlinEq = 0; sizes.mNonlinIneq = 0; % No nonlinear constraints

display = optimget(options,'Display',defaultopt,'fast');
detailedExitMsg = ~isempty(strfind(display,'detailed'));
switch display
    case {'off','none'}
        verbosity = 0;
    case {'iter','iter-detailed'}
        verbosity = 2;
    case {'final','final-detailed'}
        verbosity = 1;
    case 'testing'
        verbosity = Inf;
    otherwise
        verbosity = 1;
end
diagnostics = strcmpi(optimget(options,'Diagnostics',defaultopt,'fast'),'on');
gradflag =  strcmp(optimget(options,'Jacobian',defaultopt,'fast'),'on');

funValCheck = strcmp(optimget(options,'FunValCheck',defaultopt,'fast'),'on');

algorithm = optimget(options,'Algorithm',defaultopt,'fast');
initDamping = optimget(options,'InitDamping',defaultopt,'fast');
if ~iscell(algorithm)
    initLMparam = initDamping;
else
    initLMparam = algorithm{2}; % Initial Levenberg-Marquardt parameter
    algorithm = algorithm{1};   % Algorithm string
end

switch algorithm
    case 'trust-region-dogleg'
         algorithmflag = 2;
    case 'trust-region'
        algorithmflag = 1;
    case 'levenberg-marquardt'
        algorithmflag = 3;
    otherwise % Invalid choice of Algorithm
        error(message('optim:fsolve:InvalidAlgorithm'))
end

% Process user function
if ~isempty(FUN)  % will detect empty string, empty matrix, empty cell array
    funfcn = lsqfcnchk(FUN,'fsolve',length(varargin),funValCheck,gradflag);
else
    error(message('optim:fsolve:InvalidFUN'))
end

mtxmpy = optimget(options,'JacobMult',defaultopt,'fast');
% Check if name clash
functionNameClashCheck('JacobMult',mtxmpy,'atamult','optim:fsolve:JacobMultNameClash');

% Use internal Jacobian-multiply function if user does not provide JacobMult function 
% or options.Jacobian is off
if isempty(mtxmpy) || (~strcmpi(funfcn{1},'fungrad') && ~strcmpi(funfcn{1},'fun_then_grad'))
    mtxmpy = @atamult;
end

JAC = [];
x(:) = xstart;
switch funfcn{1}
    case 'fun'
        try
            fuser = feval(funfcn{3},x,varargin{:});
        catch userFunExcept
            optimExcept = MException('optim:fsolve:ObjectiveError', ...
               getString(message('optim:fsolve:ObjectiveError')));
            userFunExcept = addCause(userFunExcept,optimExcept);
            rethrow(userFunExcept)
        end
    case 'fungrad'
        try
            [fuser,JAC] = feval(funfcn{3},x,varargin{:});
        catch userFunExcept
            optimExcept = MException('optim:fsolve:ObjectiveError', ...
                getString(message('optim:fsolve:ObjectiveError')));
            userFunExcept = addCause(userFunExcept,optimExcept);
            rethrow(userFunExcept)
        end
    case 'fun_then_grad'
        try
            fuser = feval(funfcn{3},x,varargin{:});
        catch userFunExcept
            optimExcept = MException('optim:fsolve:ObjectiveError', ...
                getString(message('optim:fsolve:ObjectiveError')));
            userFunExcept = addCause(userFunExcept,optimExcept);
            rethrow(userFunExcept)
        end
        try
            JAC = feval(funfcn{4},x,varargin{:});
        catch userFunExcept
            optimExcept = MException('optim:fsolve:JacobianError', ...
                getString(message('optim:fsolve:JacobianError')));
            userFunExcept = addCause(userFunExcept,optimExcept);
            rethrow(userFunExcept)
        end
    otherwise
        error(message('optim:fsolve:UndefinedCalltype'))
end

% Check for non-double data typed values returned by user functions 
if ~isempty( isoptimargdbl('FSOLVE', {'F','J'}, fuser, JAC) )
    error('optim:fsolve:NonDoubleFunVal',getString(message('optimlib:commonMsgs:NonDoubleFunVal','FSOLVE')));
end

f = fuser(:);
sizes.nFun = length(f);

if gradflag
    % check size of JAC
    [Jrows, Jcols] = size(JAC);
    if isempty(options.JacobMult)
        % Not using 'JacobMult' so Jacobian must be correct size
        if Jrows ~= sizes.nFun || Jcols ~= sizes.nVar
            error(message('optim:fsolve:InvalidJacobian', sizes.nFun, sizes.nVar))
        end
    end
else
    Jrows = sizes.nFun;
    Jcols = sizes.nVar;
end

caller = 'fsolve';

% Choose what algorithm to run: determine algorithmflag and check criteria
if algorithmflag == 1 && sizes.nFun < sizes.nVar
    % trust-region algorithm and not enough equations - switch to
    % levenberg-marquardt algorithm
    warning(message('optim:fsolve:FewerFunsThanVars'))
    algorithmflag = 3;
elseif algorithmflag == 2 && sizes.nFun ~= sizes.nVar
    warning(message('optim:fsolve:NonSquareSystem'));
    algorithmflag = 3;
end

% Set up for diagnostics and derivative check
confcn = {''};
if diagnostics
    % Do diagnostics on information so far
    constflag = false; gradconstflag = false; hessflag = false; 
    non_eq = 0;non_ineq = 0;lin_eq = 0;lin_ineq = 0;
   
    % Set OUTPUT.algorithm for diagnostics
     switch algorithmflag
         case 1
             OUTPUT.algorithm = 'trust-region';
         case 2
             OUTPUT.algorithm = 'trust-region-dogleg';
         case 3
             OUTPUT.algorithm = 'levenberg-marquardt';
    end
    diagnose('fsolve',OUTPUT,gradflag,hessflag,constflag,gradconstflag,...
        xstart,non_eq,non_ineq,lin_eq,lin_ineq,LB,UB,funfcn,confcn);
end

% Read options for finitedifferences
options.FinDiffType = optimget(options,'FinDiffType',defaultopt,'fast');
options.GradObj = optimget(options,'Jacobian',defaultopt,'fast');
options.GradConstr = 'off';
DerivativeCheck = strcmp(optimget(options,'DerivativeCheck',defaultopt,'fast'),'on');
% Read in and error check option TypicalX
[typicalx,ME] = getNumericOrStringFieldValue('TypicalX','ones(numberOfVariables,1)', ...
    ones(sizes.nVar,1),'a numeric value',options,defaultopt);
if ~isempty(ME)
    throw(ME)
end
checkoptionsize('TypicalX', size(typicalx), sizes.nVar);
options.TypicalX = typicalx;
options.DiffMinChange = optimget(options,'DiffMinChange',defaultopt,'fast');
options.DiffMaxChange = optimget(options,'DiffMaxChange',defaultopt,'fast');
options = validateFinDiffRelStep(sizes.nVar,options,defaultopt);
options.UseParallel = optimget(options,'UseParallel',defaultopt,'fast');

% Create structure of flags for finitedifferences
finDiffFlags.fwdFinDiff = strcmpi(options.FinDiffType,'forward'); % Check for forward fin-diff
finDiffFlags.scaleObjConstr = false; % No scaling
finDiffFlags.chkComplexObj = false;  % Don't check whether objective function values are complex
finDiffFlags.isGrad = false;         % Compute Jacobian, not gradient
finDiffFlags.hasLBs = false(sizes.nVar,1);  % No bounds
finDiffFlags.hasUBs = false(sizes.nVar,1);
% finDiffFlags.chkFunEval will be set below according to 
% where finDiffFlags.chkFunEval will be used.

% Check derivatives
if DerivativeCheck && gradflag          % user wants to check derivatives
    lb = -Inf(sizes.nVar,1); ub = Inf(sizes.nVar,1);
    finDiffFlags.chkFunEval = false;    % don't check function values during finite differences
    validateFirstDerivatives(funfcn,confcn,x, ...
        lb,ub,options,finDiffFlags,sizes,varargin{:});
end

% Execute algorithm
if algorithmflag == 1   % trust-region
    if ~gradflag
        Jstr = optimget(options,'JacobPattern',defaultopt,'fast');
        if ischar(Jstr)
            % options.JacobPattern is the default: 'sparse(ones(jrows,jcols))'
            Jstr = sparse(ones(Jrows,Jcols));
        end
        checkoptionsize('JacobPattern', size(Jstr), Jcols, Jrows);
    else
        Jstr = [];
    end
    computeLambda = 0;
    % Set MaxFunEvals appropriately for trust-region
    defaultopt.MaxFunEvals = '100*numberOfVariables';
    % Don't check function values during finite differences
    finDiffFlags.chkFunEval = false;     
    [x,FVAL,~,JACOB,EXITFLAG,OUTPUT,msgData]=...
        snls(funfcn,x,LB,UB,verbosity,options,defaultopt,f,JAC,caller,Jstr,...
        computeLambda,mtxmpy,detailedExitMsg,optionFeedback,finDiffFlags,varargin{:});
    
    % Correct the algorithm name stored in snls
    OUTPUT.algorithm = 'trust-region';
elseif algorithmflag == 2   % trust-region dogleg    
    % Set MaxFunEvals appropriately for trust-region-dogleg
    defaultopt.MaxFunEvals = '100*numberOfVariables';
    % Check function values during finite differences
    finDiffFlags.chkFunEval = true;    
    
    [x,FVAL,JACOB,EXITFLAG,OUTPUT,msgData]=...
        trustnleqn(funfcn,x,verbosity,gradflag,options,defaultopt,f,JAC,...
        detailedExitMsg,optionFeedback,finDiffFlags,sizes,varargin{:});
elseif algorithmflag == 3   % Levenberg-Marquardt
    % Set MaxFunEvals appropriately for LM
    defaultopt.MaxFunEvals = '200*numberOfVariables';
    % Check function values during finite differences
    finDiffFlags.chkFunEval = true;    
    
    [x,FVAL,JACOB,EXITFLAG,OUTPUT,msgData] = ...
        levenbergMarquardt(funfcn,x,verbosity,options,defaultopt,f,JAC,caller, ...
        initLMparam,detailedExitMsg,optionFeedback,finDiffFlags,varargin{:});
end

Resnorm = FVAL'*FVAL;  % assumes FVAL still a vector
sqrtTolFunValue = sqrt(optimget(options,'TolFunValue',defaultopt,'fast'));
if EXITFLAG > 0 % if we think we converged:
    % Call createExitMsg with appended additional information on the closeness
    % to a root.
    if Resnorm > sqrtTolFunValue
        msgData = internalFlagForExitMessage(algorithmflag == 2,msgData,EXITFLAG);
        EXITFLAG = -2;
    end  
    %OUTPUT.message = createExitMsg(msgData{:},Resnorm,optionFeedback.TolFunValue,sqrtTolFunValue);
    OUTPUT.message = msgData;
else
    %OUTPUT.message = createExitMsg(msgData{:});
    OUTPUT.message = msgData;
end

% Reset FVAL to shape of the user-function output, fuser
FVAL = reshape(FVAL,size(fuser));
end

function msgData = internalFlagForExitMessage(isDogleg,msgData,EXITFLAG)
% internalFlagForExitMessage Create internal flag for exit message.
% Change internal exitflag to unique identifier -21, -22, -23 or -27

if isDogleg && msgData{2} == 26
    % Algorithm dogleg, stopped because of undefined Jacobian
    msgData{2} = -27;
else
    % All algorithms: change internal exitflag by negating the
    % exitflag and adding to -20.
    msgData{2} = -20 - EXITFLAG;
end
end

function [x,Fvec,JAC,EXITFLAG,OUTPUT,msgData]= trustnleqn(funfcn,x,verbosity,gradflag, ...
  options,defaultopt,Fvec,JAC,detailedExitMsg,optionFeedback,finDiffFlags,sizes,varargin)
%TRUSTNLEQN Trust-region dogleg nonlinear systems of equation solver.
%
%   TRUSTNLEQN solves a system of nonlinear equations using a dogleg trust
%   region approach.  The algorithm implemented is similar in nature
%   to the FORTRAN program HYBRD1 of J.J. More', B.S. Garbow and K.E. 
%   Hillstrom, User Guide for MINPACK 1, Argonne National Laboratory, 
%   Rept. ANL-80-74, 1980, which itself was based on the program CALFUN 
%   of M.J.D. Powell, A Fortran subroutine for solving systems of
%   nonlinear algebraic equations, Chap. 7 in P. Rabinowitz, ed.,
%   Numerical Methods for Nonlinear Algebraic Equations, Gordon and
%   Breach, New York, 1970.

%   Copyright 1990-2015 The MathWorks, Inc.
%
% NOTE: 'x' passed in and returned in matrix form.
%       'Fvec' passed in and returned in vector form.
%
% Throughout this routine 'x' and 'F' are matrices while
% 'xvec', 'xTrial', 'Fvec' and 'FTrial' are vectors. 
% This was done for compatibility with the 'fsolve.m' interface.

% Check to see if the function vector and user-supplied Jacobian at the
% initial point are Inf or NaN. If not, then terminate immediately.
% NOTE: complex values are ok in equation solving.
if any(~isfinite(Fvec)) 
    error(message('optim:trustnleqn:UsrObjUndefAtX0'));    
end
if any(~isfinite(nonzeros(JAC)))
    caller = 'fsolve';
    error('optimlib:trustnleqn:UserJacUndefAtX0', ...
        getString(message('optimlib:commonMsgs:JacUndefAtX0',caller)));
end

% Define some sizes.
xvec = x(:);         % vector representation of x
% Convert values to full to avoid unnecessary sparse operation overhead
Fvec = full(Fvec); 

% Get user-defined options.
[maxfunc,maxit,tolf,tolfunvalue,tolx,giventypx,outputfcn,plotfcns] = ...
    getOpts(sizes.nVar,options,defaultopt);

typx = options.TypicalX;
scale = giventypx;  % scaling featured only enabled when typx values provided

% For parallel finite difference (if needed) we need to send the function
% handles now to the workers. This avoids sending the function handles in
% every iteration of the solver. The output from 'setOptimFcnHandleOnWorkers'
% is a onCleanup object that will perform cleanup task on the workers.
UseParallel = optimget(options,'UseParallel',defaultopt,'fast');
cleanupObj = setOptimFcnHandleOnWorkers(UseParallel,funfcn,{''}); %#ok<NASGU>

% Handle the output function.
if isempty(outputfcn)
    haveoutputfcn = false;
else
    haveoutputfcn = true;
    xOutputfcn = x; % Last x passed to outputfcn; has the input x's shape
    % Parse OutputFcn which is needed to support cell array syntax for OutputFcn.
    outputfcn = createCellArrayOfFunctions(outputfcn,'OutputFcn');
end

% Handle the plot function.
if isempty(plotfcns)
    haveplotfcn = false;
else
    haveplotfcn = true;
    xOutputfcn = x; % Last x passed to outputfcn; has the input x's shape
    % Parse PlotFcns which is needed to support cell array syntax for PlotFcns.
    plotfcns = createCellArrayOfFunctions(plotfcns,'PlotFcns');
end

% Initialize local arrays.
d       = zeros(sizes.nVar,1);
scalMat = ones(sizes.nVar,1); 

% Initialize some trust region parameters.
Delta    = 1e0;
DeltaMax = 1e10;
eta1     = 0.05;
eta2     = 0.9;
alpha1   = 2.5;
alpha2   = 0.25;

% Other initializations.
iter = 0;
numFevals = 1;   % computed in fsolve.m
if gradflag
  numJevals = 1; % computed in fsolve.m
else
  numJevals = 0;
end 
stepAccept = true;
normd = 0.0e0;
scalemin = eps;
scalemax = 1/scalemin;
objold = 1.0e0;
obj = 0.5*(Fvec'*Fvec);  % Initial Fvec computed in fsolve.m

% Initialize the boolean indicating whether all function values are well
% defined during estimation of Jacobian.
evalOK = true;

% Define empty constraints for computeFinDiffGradAndJac
confcn = {'','','','',''};

numFDfevals = 0;
% Compute initial finite difference Jacobian, objective and gradient.
if ~gradflag
    % User is not specifying the Jacobian. Compute finite differences of
    % full Jacobian at initial point.
    JACfindiff = zeros(sizes.nFun,sizes.nVar); % pre-allocate derivative array
    % Set flags so that finite differences will validate gradient at x0
    % (only at x0, no validation of gradients a subsequent iterates).
    finDiffFlagsX0 = finDiffFlags; 
    finDiffFlagsX0.chkFunEval = true; % validate function values for grad at x0
   
    [JACfindiff,~,~,numFDfevals,evalOK] = computeFinDiffGradAndJac(x,funfcn,confcn,Fvec, ...
        [],[],JACfindiff,[],[],-Inf*ones(sizes.nVar,1),Inf*ones(sizes.nVar,1),[],options,finDiffFlagsX0,sizes,varargin{:});
    
    if ~evalOK
        caller = 'fsolve';        
        error('optim:trustnleqn:DerivUndefAtX0', ...
            getString(message('optimlib:commonMsgs:FinDiffJacUndefAtX0',caller)));
    end
    
    % Set initial Jacobian to be the finite difference approximation.
    JAC = JACfindiff;
end

% Increment the number of function evaluations with those used in the
% finite difference calls.
numFevals = numFevals + numFDfevals;

grad = JAC'*Fvec;
normgradinf = norm(grad,inf);

% Print header.
header = sprintf(['\n                                         Norm of      First-order   Trust-region\n',...
                    ' Iteration  Func-count     f(x)          step         optimality    radius']);
formatstr = ' %5.0f      %5.0f   %13.6g  %13.6g   %12.3g    %12.3g';
if verbosity > 1
  disp(header);
end

% Initialize the output function.
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,xvec,xOutputfcn,'init',iter, ...
        numFevals,Fvec,[],[],[],Delta,stepAccept,varargin{:});
    if stop
        [x,Fvec,JAC,EXITFLAG,OUTPUT] = cleanUpInterrupt(xOutputfcn,optimValues);
        msgData = {'trustnleqn',EXITFLAG,verbosity > 0,detailedExitMsg,'fsolve'};
        return;
    end
end

% Compute initial diagonal scaling matrix.
if scale
  if giventypx && ~isempty(typx) % scale based on typx values
    typx(typx==0) = 1; % replace any zero entries with ones
    scalMat = 1./abs(typx);
  else         % scale based on norm of the Jacobian (not currently active)  
    scalMat = getscalMat(sizes.nVar,JAC,scalemin,scalemax);
  end
end

% Display initial iteration information.
formatstr0 = ' %5.0f      %5.0f   %13.6g                  %12.3g    %12.3g';
% obj is 0.5*F'*F but want to display F'*F
iterOutput0 = sprintf(formatstr0,iter,numFevals,2*obj,normgradinf,Delta);
if verbosity > 1
   disp(iterOutput0);
end
% OutputFcn call.
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,xvec,xOutputfcn,'iter',iter, ...
        numFevals,Fvec,normd,grad,normgradinf,Delta,stepAccept,varargin{:});
    if stop
        [x,Fvec,JAC,EXITFLAG,OUTPUT] = cleanUpInterrupt(xOutputfcn,optimValues);
        msgData = {'trustnleqn',EXITFLAG,verbosity > 0,detailedExitMsg,'fsolve'};
        return;
    end
end

% Test convergence at initial point.
[done,EXITFLAG,msgData] = testStop(normgradinf,tolf,tolfunvalue,tolx,...
     stepAccept,iter,evalOK,maxit,numFevals,maxfunc,Delta,normd,...
     obj,objold,d,xvec,detailedExitMsg,optionFeedback,verbosity);

% Create the fault tolerance structure
faultTolStruct = createFaultTolStruct(false);
newFaultTolStruct = faultTolStruct;

% Beginning of main iteration loop.
while ~done
  iter = iter + 1;

  % Compute step, d, using dogleg approach.
  [d,quadObj,normd,normdscal] = ...
       dogleg(sizes.nVar,Fvec,JAC,grad,Delta,scalMat,varargin);

  % Compute the model reduction given by d (pred).
  pred = -quadObj;

  % Compute the trial point, xTrial.
  xTrial = xvec + d;

  % Evaluate nonlinear equations and objective at trial point.
  switch funfcn{1}
  case 'fun'
    F = feval(funfcn{3},reshape(xTrial,sizes.xShape),varargin{:});
  case 'fungrad'
    [F,JACTrial] = feval(funfcn{3},reshape(xTrial,sizes.xShape),varargin{:});
    numJevals = numJevals + 1;
  case 'fun_then_grad'
    F = feval(funfcn{3},reshape(xTrial,sizes.xShape),varargin{:}); 
  otherwise
    error(message('optim:trustnleqn:UndefinedCalltype'))
  end  
  numFevals = numFevals + 1;
  FTrial = full(F(:)); % make FTrial a vector, convert to full
  objTrial = 0.5*(FTrial'*FTrial); 

  % Compute the actual reduction given by xTrial (ared).
  ared = obj - objTrial;

  % Compute ratio = ared/pred.
  if pred <= 0 % reject step
    ratio = 0;
  else
    ratio = ared/pred;
  end
  
  % Update fault tolerance structure.
  faultTolStruct = updateFaultTolStruct(faultTolStruct, objTrial, ...
      verbosity > 1);
  
  if haveoutputfcn % Call output functions (we don't call plot functions with 'interrupt' flag)
      [~, ~, stop] = callOutputAndPlotFcns(outputfcn,{},xvec,xOutputfcn,'interrupt',iter, ...
          numFevals,Fvec,normd,grad,normgradinf,Delta,stepAccept,varargin{:});
      if stop  % Stop per user request.
          [x,Fvec,JAC,EXITFLAG,OUTPUT] = cleanUpInterrupt(xOutputfcn,optimValues);
          msgData = {'trustnleqn',EXITFLAG,verbosity > 0,detailedExitMsg,'fsolve'};
          return;
      end
  end
  
  if ratio > eta1 && faultTolStruct.currTrialWellDefined % accept step.

    xvec = xTrial; Fvec = FTrial; objold = obj; obj = objTrial;
    x(:) = xvec; % update matrix representation
    % Compute JAC at new point. (already computed with F if 'fungrad')

    % Compute finite difference Jacobian if needed.
    if ~gradflag

        [JACfindiff,~,~,numFDfevals,evalOK] = computeFinDiffGradAndJac(x,funfcn,confcn,Fvec, ...
            [],[],JACfindiff,[],[],-Inf*ones(sizes.nVar,1),Inf*ones(sizes.nVar,1),[],options,finDiffFlags,sizes,varargin{:});
        
        numFevals = numFevals + numFDfevals;
    end

    switch funfcn{1}
        case 'fun'
            JAC = JACfindiff;
        case 'fungrad'
            JAC = JACTrial;
        case 'fun_then_grad'
            JAC = feval(funfcn{4},x,varargin{:});
            numJevals = numJevals + 1;
        otherwise
            error(message('optim:trustnleqn:UndefinedCalltype'))
    end
      
    grad = JAC'*Fvec;
    normgradinf = norm(grad,inf);

    % Update internal diagonal scaling matrix (dynamic scaling).
    if scale && ~giventypx
      scalMat = getscalMat(sizes.nVar,JAC,scalemin,scalemax);
    end

    stepAccept = true;
  else % reject step.
    stepAccept = false;
  end 

  % Print iteration statistics.
  if verbosity > 1
      if faultTolStruct.undefObj
          fprintf(getString(message('optimlib:commonMsgs:ObjInfNaNComplex', ...
              faultTolStruct.undefValue)));
      end
      % obj is 0.5*F'*F but want to display F'*F
      iterOutput = sprintf(formatstr,iter,numFevals,2*obj,normd,normgradinf,Delta);
      disp(iterOutput);
  end
  % OutputFcn call.
  if haveoutputfcn || haveplotfcn
      [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,xvec,xOutputfcn,'iter',iter, ...
        numFevals,Fvec,normd,grad,normgradinf,Delta,stepAccept,varargin{:});
      if stop
          [x,Fvec,JAC,EXITFLAG,OUTPUT] = cleanUpInterrupt(xOutputfcn,optimValues);
          msgData = {'trustnleqn',EXITFLAG,verbosity > 0,detailedExitMsg,'fsolve'};
          return;
      end
  end

  % Update trust region radius.
  Delta = updateDelta(Delta,ratio,normdscal,eta1,eta2,...
                      alpha1,alpha2,DeltaMax,faultTolStruct);

  % Check for termination.
  [done,EXITFLAG,msgData] = testStop(normgradinf,tolf,tolfunvalue,tolx,...
       stepAccept,iter,evalOK,maxit,numFevals,maxfunc,Delta,normd,...
       obj,objold,d,xvec,detailedExitMsg,optionFeedback,verbosity);
   
  % As the iteration has been completed, the fault tolerance
  % structure needs to be reset.
  faultTolStruct = newFaultTolStruct;
   
end

if haveoutputfcn || haveplotfcn
    callOutputAndPlotFcns(outputfcn,plotfcns,xvec,xOutputfcn,'done',iter, ...
        numFevals,Fvec,normd,grad,normgradinf,Delta,stepAccept,varargin{:});
    % Optimization done, so ignore "stop"
end


% Optimization is finished.

% Assign output statistics.
OUTPUT.iterations = iter;
OUTPUT.funcCount = numFevals;
OUTPUT.algorithm = 'trust-region-dogleg';
OUTPUT.firstorderopt = normgradinf;
end
% TRUSTNLEQN finished

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [maxfunc,maxit,tolf,tolfunvalue,tolx,giventypx,outputfcn,plotfcns] = ...
                getOpts(nvar,options,defaultopt)
%getOpts gets the user-defined options for TRUSTNLEQN.

% Both Medium and Large-Scale options.
maxfunc = optimget(options,'MaxFunEvals',defaultopt,'fast');
if ischar(maxfunc)
  if isequal(lower(maxfunc),'100*numberofvariables')
    maxfunc = 100*nvar;
  else
    error(message('optim:trustnleqn:InvalidMaxFunEvals'))
  end
end
maxit = optimget(options,'MaxIter',defaultopt,'fast');
tolfunvalue = optimget(options,'TolFunValue',defaultopt,'fast');
tolf = optimget(options,'TolFun',defaultopt,'fast');
tolx = optimget(options,'TolX',defaultopt,'fast');
outputfcn = optimget(options,'OutputFcn',defaultopt,'fast');
plotfcns = optimget(options,'PlotFcns',defaultopt,'fast');

% Check if TypicalX is the default (ones) 
giventypx = any(options.TypicalX ~= 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [done,EXITFLAG,msgData] = testStop(normgradinf,tolf,tolfunvalue,...
    tolx,stepAccept,iter,evalOK,maxit,numFevals,maxfunc,Delta,normd,...
    obj,objold,d,xvec,detailedExitMsg,optionFeedback,verbosity)
%testStop checks the termination criteria for TRUSTNLEQN.

done = false;
EXITFLAG = 0;
msgData = {};

% Check termination criteria.
if ~evalOK
    EXITFLAG = 2; 
    msgFlag = 26;
    msgData = {'trustnleqn',msgFlag,verbosity > 0,detailedExitMsg,'fsolve', ...
        [], [], [], [], [], []};
    done = true;
elseif stepAccept && normgradinf < tolf
  done = true;
  EXITFLAG = 1;
  if iter == 0
      msgFlag = 100;
  else
      msgFlag = EXITFLAG;
  end
  % Setup input parameters for createExitMsg with msgFlag = 100 if x0 is
  % optimal, otherwise msgFlag = 1
  msgData = {'trustnleqn',msgFlag,verbosity > 0,detailedExitMsg,'fsolve', ...
      normgradinf,optionFeedback.TolFun,tolf,2*obj,optionFeedback.TolFun,sqrt(tolf)};
elseif iter > 1 && max(abs(d)./(abs(xvec)+1)) < max(tolx^2,eps)
   % Assign msgFlag, a unique internal exitflag, to 2 or -22 for this
   % stopping test depending on whether the result appears to be a root or
   % not.
   if 2*obj < sqrt(tolf) % fval'*fval < sqrt(tolf)
      EXITFLAG = 2; msgFlag = 2;
      dispMsg = verbosity > 0;
   else
      EXITFLAG = -2; msgFlag = -22;
      dispMsg = verbosity > 0;
   end
   % Setup input parameters for createExitMsg
   msgData = {'trustnleqn',msgFlag,dispMsg,detailedExitMsg,'fsolve', ...
       max(abs(d)./(abs(xvec)+1)),optionFeedback.TolX,max(tolx^2,eps), ...
       2*obj,optionFeedback.TolFun,sqrt(tolf)};
   done = true;
elseif iter > 1 && stepAccept && normd < 0.9*Delta ...
                && abs(objold-obj) < max(tolfunvalue^2,eps)*(1+abs(objold))
  % Assign msgFlag, a unique internal exitflag, to 3 or -23 for this
  % stopping test depending on whether the result appears to be a root or
  % not.
  if 2*obj < sqrt(tolfunvalue) % fval'*fval < sqrt(tolf)
     EXITFLAG = 3; msgFlag = 3;
     dispMsg = verbosity > 0;
  else
     EXITFLAG = -2; msgFlag = -23;
     dispMsg = verbosity > 0;
  end
  % Setup input parameters for createExitMsg
  msgData = {'trustnleqn',msgFlag,dispMsg,detailedExitMsg,'fsolve', ...
      abs(objold-obj)./(abs(objold)+1),optionFeedback.TolFunValue,...
      max(tolfunvalue^2,eps),2*obj,optionFeedback.TolFunValue,sqrt(tolfunvalue)};
  done = true;
elseif Delta < 2*eps
  EXITFLAG = -3;
  msgData = {'trustnleqn',EXITFLAG,verbosity > 0,detailedExitMsg,'fsolve', ...
      Delta,'',2*eps};
  done = true;
elseif iter >= maxit
  EXITFLAG = 0;
  msgData = {'trustnleqn',10,verbosity > 0,detailedExitMsg,'fsolve', ...
      [],optionFeedback.MaxIter,maxit};
  done = true;
elseif numFevals >= maxfunc
  EXITFLAG = 0;
  msgData = {'trustnleqn',EXITFLAG,verbosity > 0,detailedExitMsg,'fsolve', ...
      [],optionFeedback.MaxFunEvals,maxfunc};
  done = true;
end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Delta = updateDelta(Delta,ratio,normdscal,eta1,eta2,...
                             alpha1,alpha2,DeltaMax,faultTolStruct)
%updateDelta updates the trust region radius in TRUSTNLEQN.
%
%   updateDelta updates the trust region radius based on the value of
%   ratio and the norm of the scaled step.

if ~faultTolStruct.currTrialWellDefined
    % Shrink the trust region if any element of the function vector
    % at the new step is not defined (i.e. it's inf/NaN/complex). The
    % update for Delta in this case follows that in snls.
    Delta = min(normdscal/20,Delta/20);
elseif ratio < eta1
    Delta = alpha2*normdscal;
elseif ratio >= eta2
    Delta = max(Delta,alpha1*normdscal);
end
Delta = min(Delta,DeltaMax);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scalMat = getscalMat(nvar,JAC,scalemin,scalemax)
%getscalMat computes the scaling matrix in TRUSTNLEQN.
%
%   getscalMat computes the scaling matrix based on the norms 
%   of the columns of the Jacobian.

scalMat = ones(nvar,1);
for i=1:nvar
  scalMat(i,1) = norm(JAC(:,i));
end
scalMat(scalMat<scalemin) = scalemin;  % replace small entries
scalMat(scalMat>scalemax) = scalemax;  % replace large entries
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
function [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,xvec,xOutputfcn,state,iter,numFevals, ...
    Fvec,normd,grad,normgradinf,Delta,stepAccept,varargin)
% CALLOUTPUTANDPLOTFCNS assigns values to the struct OptimValues and then calls the
% outputfcn/plotfcns.  
%
% state - can have the values 'init','iter','interrupt', or 'done'. 
%
% For the 'done' state we do not check the value of 'stop' because the
% optimization is already done.

optimValues.iteration = iter;
optimValues.funccount = numFevals;
optimValues.fval = Fvec;
optimValues.stepsize = normd; 
optimValues.gradient = grad; 
optimValues.firstorderopt = normgradinf;
optimValues.trustregionradius = Delta;
optimValues.stepaccept = stepAccept;

xOutputfcn(:) = xvec;  % Set xvec to have user expected size
stop = false;
% Call output functions
if ~isempty(outputfcn)
    switch state
        case {'iter','init','interrupt'}
            stop = callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error(message('optim:trustnleqn:UnknownStateInCALLOUTPUTANDPLOTFCNS'))
    end
end
% Call plot functions
if ~isempty(plotfcns)
    switch state
        case {'iter','init'}
            stop = callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error(message('optim:trustnleqn:UnknownStateInCALLOUTPUTANDPLOTFCNS'))
    end
end
end
%--------------------------------------------------------------------------
function [x,Fvec,JAC,EXITFLAG,OUTPUT] = cleanUpInterrupt(xOutputfcn,optimValues)
% CLEANUPINTERRUPT sets the outputs arguments to be the values at the last call
% of the outputfcn during an 'iter' call (when these values were last known to
% be consistent). 

% Call plot function driver to finalize the plot function figure window. If
% no plot functions have been specified or the plot function figure no
% longer exists, this call just returns.
callAllOptimPlotFcns('cleanuponstopsignal');

x = xOutputfcn; 
Fvec = optimValues.fval;
EXITFLAG = -1; 
OUTPUT.iterations = optimValues.iteration;
OUTPUT.funcCount = optimValues.funccount;
OUTPUT.algorithm = 'trust-region-dogleg';
OUTPUT.firstorderopt = optimValues.firstorderopt; 
JAC = []; % May be in an inconsistent state
end

function [step,quadObj,normStep,normStepScal] = ...
    dogleg(nvar,F,JAC,grad,Delta,scalMat,varargin)
%DOGLEG approximately solves trust region subproblem via a dogleg approach.
%
%   DOGLEG finds an approximate solution d to the problem:
%
%     min_d      f + g'd + 0.5*d'Bd
%
%     subject to ||Dd|| <= Delta
%
%   where g is the gradient of f, B is a Hessian approximation of f, D is a
%   diagonal scaling matrix and Delta is a given trust region radius.

%   Copyright 1990-2012 The MathWorks, Inc.

% NOTE: The scaling matrix D above is called scalMat in this routine.

% Compute scaled gradient and other scaled terms.
gradscal = grad./scalMat;
gradscal2 = gradscal./scalMat;
normgradscal = norm(gradscal);

if normgradscal >= eps
    % First compute the Cauchy step (in scaled space).
    dCauchy = -(Delta/normgradscal)*gradscal;
    JACvec = JAC*gradscal2;
    denom = Delta*(JACvec'*JACvec);
    tauterm = normgradscal^3/denom;
    tauC = min(1,tauterm);
    dCauchy = tauC*dCauchy;
    
    % Compute quadratic objective at Cauchy point.
    JACvec = JAC*(dCauchy./scalMat); 
    objCauchy = gradscal'*dCauchy + 0.5*(JACvec'*JACvec);   
    normdCauchy = min(norm(dCauchy),Delta);
else
    % Set Cauchy step to zero step and continue.
    objCauchy = 0;
    normdCauchy = 0;
    dCauchy = zeros(nvar,1);
end

if Delta - normdCauchy < eps;
    % Take the Cauchy step if it is at the boundary of the trust region.
    step = dCauchy; quadObj = objCauchy;
else
    % Compute the Gauss-Newton step (in scaled space).
    % Disable the warnings about conditioning for singular and
    % nearly singular matrices
    warningstate1 = warning('off', 'MATLAB:nearlySingularMatrix');
    warningstate2 = warning('off', 'MATLAB:singularMatrix');
    dNewton = -JAC\F;
    % Restore the warning states to their original settings
    warning(warningstate1)
    warning(warningstate2)
    
    dNewton = dNewton.*scalMat;     % scale the step
    
    if any(~isfinite(dNewton))
        % Take the Cauchy step if the Gauss-Newton step gives bad values.
        step = dCauchy; quadObj = objCauchy;
    else 
        normdNewt = norm(dNewton);
        if normdNewt <= Delta
            % Use the Newton direction as the trial step
            step = dNewton;
        else
            % Find the intersect point along dogleg path.
            Delta2 = Delta^2;
            normdCauchy2 = min(normdCauchy^2,Delta2);
            normdNewt2 = normdNewt^2;
            dCdN = dCauchy'*dNewton;
            dCdNdist2 = max((normdCauchy2+normdNewt2-2*dCdN),0);
            
            if dCdNdist2 == 0
                tauI = 0;
            else
                % Stable method for solving 1-D quadratic
                a = 0.5*dCdNdist2;
                b = dCdN - normdCauchy2;
                c = 0.5*(normdCauchy2 - Delta2);
                q = -0.5*(b + sign(b)*sqrt(b^2 - 4*a*c));
                if b > 0
                    tauI = c/q; 
                else
                    tauI = q/a;
                end

                % Make sure we take a finite step,
                % (check for poorly scaled/infinite directions).
                if ~isfinite(tauI)
                    tauI = 0; % Take Cauchy step
                end
            end
            step = dCauchy + tauI*(dNewton-dCauchy);
        end
        % Compute quadratic objective at trial point.
        JACvec = JAC*(step./scalMat);
        quadObj = gradscal'*step + 0.5*(JACvec'*JACvec);
        
        % Compare Cauchy step and trial step (Newton or intersection)
        if objCauchy < quadObj
            step = dCauchy; quadObj = objCauchy;
        end
    end
end

% The step computed was the scaled step.  Unscale it.
normStepScal = norm(step);
step = step./scalMat;
normStep = norm(step);
end
