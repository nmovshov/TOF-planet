function tof = double_polytrope(N, x, zstrat)
%DOUBLE_POLYTROPE Model planet approximated by two polytropes.
%    DOUBLE_POLYTROPE(N, x) returns an N-point TOFPlanet object with two
%    barotropes.Polytrope eos objects. First polytrope defined by constant x(1)
%    and index x(2) is applied to levels 1:tind-1. Second polytrope defined by
%    constant x(3) and index x(4) is applied to levels tind:N. Transition is at
%    level tind, the level with (si/s0) nearest to x(5). The default level spacing
%    is one of equal radius increments between s/s0=1 and s/s0=1/N, except that
%    the transition level is snapped to the given value, x(5).
%
%    DOUBLE_POLYTROPE(N, x, zstrat) lets you specify the zvec distribution
%    strategy. Pass a handle to a function that takes a single scalar integer
%    (number of layers) and returns a vector of that length with values in the
%    interval (0, 1], for normalized mean level radii. For example, to set levels
%    with equally spaced radii use zstrat=@(n)linspace(1,1/n,n). Note that the
%    final level radii may be slightly different due to placement of the
%    transition radius. A collection of pre-made distributions is available in
%    package +zvecs.

try
    narginchk(2,3)
    if ((nargin < 3) || isempty(zstrat)), zstrat = @zvecs.best; end
    validateattributes(N, {'numeric'}, {'positive', 'integer'}, '', 'N', 1)
    validateattributes(x, {'numeric'}, {'vector', 'numel', 5, 'nonnegative'}, 2)
    validateattributes(zstrat, {'function_handle'}, {}, '', 'zstrat', 3)
    assert(x(5)>0 && x(5)<1, 'Transition (normalized) radius must be in (0,1).')
catch ME
    if nargout == 0
        help('tofmodels.double_polytrope')
    end
    rethrow(ME)
end

tof = TOFPlanet();

svec = zstrat(N);
assert(isnumeric(svec) && isvector(svec) && (numel(svec) == N),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')
assert(all(svec > 0) && all(svec <= 1),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')

[~, tind] = min(abs(svec-x(5)));
assert(tind > 2,...
    'Transition too close to surface; first polytrope has zero layers.')
svec(tind) = x(5);
tof.si = svec;
tof.rhoi = ones(N,1);

eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.Polytrope(x(3), x(4));
tof.eos = [repmat(eos1, tind - 1, 1);...
           repmat(eos2, N - tind + 1, 1)];
end
