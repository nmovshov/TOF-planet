function options = tofset(varargin)
%TOFSET Create options structure used by TOFPlanet class methods.
%   OPTIONS = TOFSET('NAME1',VALUE1,'NAME2',VALUE2,...) creates an options
%   structure OPTIONS in which the named properties have the specified values.
%   Any unspecified properties have default values. Case is ignored for property
%   names.
%
%   TOFSET with no input or output arguments displays all property names and their
%   possible values.
%
%KNOWN PROPERTIES
%
%dJtol - Convergence tolerance for gravity coefficients [ positive real {1e-6} ]
%drhotol - Convergence tolerance for density adjustment [ positive real {1e-6} ]
%MaxIterBar - Number of iterations allowed for relaxation to barotrope [ positive integer {60} ]
%MaxIterHE - Number of iterations allowed for relaxation to equilibrium shape [ positive integer {60} ]
%splineskip - Solve figure functions on every nth level and spline the rest [ nonnegative integer {0} ]
%verbosity - Level of runtime messages [0 {1} 2 3 4]
%debug - Debug mode flag (use preals instead of doubles) [ true | {false} ]

% If no arguments print usage and return.
if (nargin == 0) && (nargout == 0)
    print_usage()
    return
end

% Define name-value pairs.
p = inputParser;
p.FunctionName = mfilename;

p.addParameter('dJtol',1e-6,@isposscalar)
p.addParameter('drhotol',1e-6,@isposscalar)
p.addParameter('MaxIterBar',60,@isposintscalar)
p.addParameter('MaxIterHE',60,@isposintscalar)
p.addParameter('verbosity',1,@isnonnegintscalar)
p.addParameter('debug',false,@islogicalscalar)
p.addParameter('masmeth','layerz');
p.addParameter('splineskip',0,@isnonnegintscalar)

% Parse name-value pairs and return.
p.parse(varargin{:})
options = p.Results;

end

function isposscalar(x)
validateattributes(x,{'numeric'},{'positive','scalar'})
end

function isposintscalar(x)
validateattributes(x,{'numeric'},{'positive','integer','scalar'})
end

function isnonnegintscalar(x)
validateattributes(x,{'numeric'},{'nonnegative','integer','scalar'})
end

function islogicalscalar(x)
validateattributes(x,{'logical'},{'scalar'})
end

function print_usage()
help(mfilename)
end
