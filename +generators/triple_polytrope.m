function tof = triple_polytrope(N, x, zstrat)
%TRIPLE_POLYTROPE Model planet approximated by three polytropes.
%    TRIPLE_POLYTROPE(N, x) returns an N-point TOFPlanet object with three
%    barotropes.Polytrope eos objects. First polytrope defined by constant x(1)
%    and index x(2) is applied to levels 1:tind-1. Second polytrope defined by
%    constant x(3) and index x(4) is applied to levels tind:cind-1. Third
%    polytrope defined by constant x(5) and index x(6) is applied to levels
%    cind:N. Transition from first to second polytrope is at level tind, the level
%    with (si/s0) nearest x(7). Transition from second to third polytrope is at
%    level cind, the level with (si/s0) nearest x(8). The default level spacing is
%    one of equal radius increments between s/s0=1 and s/s0=1/N, except that the
%    transition levels are snapped to the given values, x(7) and x(8).
%
%    TRIPLE_POLYTROPE(N, x, zstrat) lets you specify the zvec distribution
%    strategy. Pass a handle to a function that takes a single scalar integer
%    (number of layers) and returns a vector of that length with values in the
%    interval (0, 1], for normalized mean level radii. For example, to set levels
%    with equally spaced radii use zstrat=@(n)linspace(1,1/n,n). Note that the
%    final level radii may be slightly different due to placement of transition
%    radii. A collection of pre-made distributions is available in package
%    +zvecs.

if nargin == 0, help('generators.triple_polytrope'), return, end
narginchk(2,3)
if ((nargin < 3) || isempty(zstrat)), zstrat = @zvecs.best; end
validateattributes(N, {'numeric'}, {'positive', 'integer'}, '', 'N', 1)
validateattributes(x, {'numeric'}, {'vector', 'numel', 8, 'nonnegative'}, 2)
validateattributes(zstrat, {'function_handle'}, {}, '', 'zstrat', 3)
assert(x(7) > 0 && x(7) < 1, 'First transition radius must be in (0,1).')
assert(x(8) > 0 && x(8) < 1, 'Second transition radius must be in (0,1).')
assert(x(8) <= x(7), 'Second transition must come before first transition.')

tof = TOFPlanet();

svec = zstrat(N);
assert(isnumeric(svec) && isvector(svec) && (numel(svec) == N),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')
assert(all(svec > 0) && all(svec <= 1),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')

[~, tind] = min(abs(svec-x(7)));
[~, cind] = min(abs(svec-x(8)));
svec(tind) = x(7);
svec(cind) = x(8);

tof.si = svec;
tof.rhoi = ones(N,1);

eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.Polytrope(x(3), x(4));
eos3 = barotropes.Polytrope(x(5), x(6));
tof.eos = [repmat(eos1, tind - 1, 1);...
           repmat(eos2, cind - tind, 1);...
           repmat(eos3, N - cind + 1, 1)];
end
