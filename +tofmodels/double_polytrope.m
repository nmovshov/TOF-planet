function tof = double_polytrope(N, x)
%DOUBLE_POLYTROPE Model planet approximated by two polytropes.
%    DOUBLE_POLYTROPE(N, x) returns an N-point TOFPlanet object with two
%    barotropes.Polytrope eos objects. First polytrope defined by constant x(1)
%    and index x(2) is applied to levels 1:tind-1. Second polytrope defined by
%    constant x(3) and index x(4) is applied to levels tind:N. Transition is at
%    level tind, the level with (si/s0) nearest to x(5). The default level spacing
%    is one of equal radius increments between s/s0=1 and s/s0=1/N, except that
%    the transition level is snapped to the given value, x(5).

try
    narginchk(2,2)
    validateattributes(N, {'numeric'}, {'positive', 'integer'}, '', 'N', 1)
    validateattributes(x, {'numeric'}, {'vector', 'numel', 5, 'nonnegative'}, 2)
    assert(x(5)>0 && x(5)<1, 'Transition (normalized) radius must be in (0,1).')
catch ME
    help('tofmodels.double_polytrope')
    rethrow(ME)
end

tof = TOFPlanet();

svec = linspace(1, 1/N, N);
[~, tind] = min(abs(svec-x(5)));
assert(tind > 2,...
    'First transition too close to surface; first polytrope has zero layers.')
svec(tind) = x(5);
tof.si = svec;
tof.rhoi = ones(N,1);

eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.Polytrope(x(3), x(4));
tof.eos = [repmat(eos1, tind - 1, 1);...
           repmat(eos2, N - tind + 1, 1)];
end
