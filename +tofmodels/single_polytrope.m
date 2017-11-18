function tof = single_polytrope(N, x, zstrat)
%SINGLE_POLYTROPE The simplest toy model planet.
%    SINGLE_POLYTROPE(N, x) returns an N-point TOFPlanet object with a
%    barotropes.Polytrope with constant x(1) and index x(2). The default level
%    spacing is one of equal radius increments between s/s0=1 and s/s0=1/N.
%
%    SINGLE_POLYTROPE(N, x, zstrat) lets you specify the zvec distribution
%    strategy. Pass a handle to a function that takes a single scalar integer
%    (number of layers) and returns a vector of that length with values in the
%    interval (0, 1], for normalized mean level radii. For example, to set levels
%    with equally spaced radii use zstrat=@(n)linspace(1,1/n,n). A collection of
%    pre-made distributions is available in package +zvecs.

try
    narginchk(2,3)
    if ((nargin < 3) || isempty(zstrat)), zstrat = @zvecs.best; end
    validateattributes(N, {'numeric'}, {'positive', 'integer'}, '', 'N', 1)
    validateattributes(x, {'numeric'}, {'vector', 'numel', 2, 'nonnegative'}, 2)
    validateattributes(zstrat, {'function_handle'}, {}, '', 'zstrat', 3)
catch ME
    help('tofmodels.single_polytrope')
    rethrow(ME)
end

tof = TOFPlanet();

svec = zstrat(N);
assert(isnumeric(svec) && isvector(svec) && (numel(svec) == N),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')
assert(all(svec > 0) && all(svec <= 1),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')

tof.si = svec;
tof.rhoi = ones(N,1);

eos1 = barotropes.Polytrope(x(1), x(2));
tof.eos = eos1;

end
