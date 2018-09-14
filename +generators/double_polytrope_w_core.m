function tof = double_polytrope_w_core(N, x, zstrat)
%DOUBLE_POLYTROPE_W_CORE Model planet with two polytropes and a solid core.
%    DOUBLE_POLYTROPE_W_CORE(N, x) returns an N-point TOFPlanet object with three
%    barotropes eos objects. The first is a polytrope defined by constant x(1) and
%    index x(2) and applied to levels 1:tind-1. Second is a polytrope defined by
%    constant x(3) and index x(4) and applied to levels tind:cind-1. Third is a
%    constant density solid defined by x(5) and applied to levels cind:N.
%    Transition from first to second polytrope is at level tind, the level with
%    (si/s0) nearest x(6). Transition from second polytrope to solid is at level
%    cind, the level with (si/s0) nearest x(7).
%
%    The default level spacing is one of equal radius increments between s/s0=1
%    and s/s0=1/N, except that the transition levels are snapped to the given
%    values, x(6) and x(7).

if nargin == 0, help('tofmodels.double_polytrope_w_core'), return, end
narginchk(2,3)
if ((nargin < 3) || isempty(zstrat)), zstrat = @zvecs.best; end
validateattributes(N, {'numeric'}, {'positive', 'integer'}, '', 'N', 1)
validateattributes(x, {'numeric'}, {'vector', 'numel', 7, 'nonnegative'}, 2)
validateattributes(zstrat, {'function_handle'}, {}, '', 'zstrat', 3)
assert(x(6) > 0 && x(6) < 1, 'Envelope transition radius must be in (0,1).')
assert(x(7) > 0 && x(7) < 1, 'Core radius must be in (0,1).')
assert(x(7) <= x(6), 'Core radius must be inside lower envelope.')

tof = TOFPlanet();

svec = zstrat(N);
assert(isnumeric(svec) && isvector(svec) && (numel(svec) == N),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')
assert(all(svec > 0) && all(svec <= 1),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')

[~, tind] = min(abs(svec-x(6)));
[~, cind] = min(abs(svec-x(7)));
svec(tind) = x(6);
svec(cind) = x(7);

tof.si = svec;
tof.rhoi = ones(N,1);

eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.Polytrope(x(3), x(4));
eos3 = barotropes.ConstDensity(x(5));
tof.eos = [repmat(eos1, tind - 1, 1);...
           repmat(eos2, cind - tind, 1);...
           repmat(eos3, N - cind + 1, 1)];
end
