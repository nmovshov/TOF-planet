function options = tofset(varargin)
%TOFSET Create options structure used by TOFPlanet class methods.
%   OPTIONS = TOFSET('NAME1',VALUE1,'NAME2',VALUE2,...) creates an options
%   structure OPTIONS in which the named properties have the specified values.
%   Any unspecified properties have default values. Case is ignored for property
%   names and unique partials are allowed.
%
%   TOFSET with no input or output arguments displays all property names and
%   their possible values.
%
%KNOWN PROPERTIES
%
%toforder - Theory of figures expansion order [ {4} 7 ]
%dJtol - Convergence tolerance for gravity coefficients [ positive real {1e-6} ]
%drhotol - Convergence tolerance for density adjustment [ positive real {1e-6} ]
%MaxIterBar - Number of iterations allowed for relaxation to barotrope [ positive integer {60} ]
%MaxIterHE - Number of iterations allowed for relaxation to equilibrium shape [ positive integer {60} ]
%xlevels - Solve figure functions on xlevels and spline the rest [ integer scalar or vector (-1 to disable) {-1} ]
%verbosity - Level of runtime messages [0 {1} 2 3 4]

% If no arguments print usage and return.
if (nargin == 0) && (nargout == 0)
    print_usage()
    return
end

% Define name-value pairs.
p = inputParser;
p.FunctionName = mfilename;

p.addParameter('toforder', 4, @(n)(n==4)||(n==7))
p.addParameter('dJtol',1e-6,@isposscalar)
p.addParameter('drhotol',1e-6,@isposscalar)
p.addParameter('MaxIterBar',60,@isposintscalar)
p.addParameter('MaxIterHE',60,@isposintscalar)
p.addParameter('xlevels',-1,@isintscalar)
p.addParameter('verbosity',1,@isnonnegintscalar)

% undocumented or obsolete options
p.addParameter('masmeth','trapz'); % undocumented mass integral method
p.addParameter('prsmeth','trapz'); % undocumented pressure integral method
p.addParameter('moimeth','midlayerz'); % undocumented moi integral method

% Parse name-value pairs and return.
p.parse(varargin{:})
options = p.Results;

end

function isposscalar(x)
validateattributes(x,{'numeric'},{'positive','scalar'})
end

function isintscalar(x)
validateattributes(x,{'numeric'},{'integer','scalar'})
end

function islogicalscalar(x) %#ok<*DEFNU>
validateattributes(x,{'logical'},{'scalar'})
end

function isposintscalar(x)
validateattributes(x,{'numeric'},{'positive','integer','scalar'})
end

function isnonnegintscalar(x)
validateattributes(x,{'numeric'},{'nonnegative','integer','scalar'})
end

function print_usage()
help(mfilename)
end
