function options = tofset(varargin)
%TOFSET Create options structure used by TOFPlanet class methods.
%   OPTIONS = TOFSET('NAME1',VALUE1,'NAME2',VALUE2,...) creates an options
%   structure OPTIONS in which the named properties have the specified values.
%   Any unspecified properties have default values. Case is ignored for property
%   names.
%
%   TOFSET with no input arguments displays all property names and their
%   possible values.
%
%KNOWN PROPERTIES
%
%dJtol - Convergence tolerance for hydrostatic equilibrium [ positive real {1e-8} ]
%dBtol - Convergence tolerance for barotrope adjustment [ positive real {1e-6} ]
%MaxIter - Maximum number of iterations allowed for relaxation [ positive integer {40} ]
%verbosity - Level of runtime messages [0 {1} 2 3 4]

% If no arguments print usage and return.
if (nargin == 0) && (nargout == 0)
    print_usage()
    return
end

% Define name-value pairs.
p = inputParser;
p.FunctionName = mfilename;

p.addParameter('dJtol',1e-8,@isposscalar)
p.addParameter('dBtol',1e-6,@isposscalar)
p.addParameter('MaxIter',40,@isposintscalar)
p.addParameter('verbosity',1,@isnonnegintscalar)

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

function TF = isvalidemail(x)
if isempty(x), TF = true; return, end
validateattributes(x,{'char'},{'row'})
validemail='[a-z_.1-9]+@[a-z_.1-9]+\.(com|net|edu)';
imatch=regexp(x,validemail);
if isempty(imatch) || ~isscalar(imatch) || imatch > 1
    error('Not a valid email address')
end
end

function print_usage()
help(mfilename)
end
