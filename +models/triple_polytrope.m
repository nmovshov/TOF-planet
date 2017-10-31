function tof = triple_polytrope(N, x)
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

try
    narginchk(2,2)
    validateattributes(N, {'numeric'}, {'positive', 'integer'}, '', 'N', 1)
    validateattributes(x, {'numeric'}, {'vector', 'numel', 8, 'nonnegative'}, 2)
    assert(x(7) > 0 && x(7) < 1, 'First transition radius must be in (0,1).')
    assert(x(8) > 0 && x(8) < 1, 'Second transition radius must be in (0,1).')
    assert(x(8) < x(7), 'Second transition must come before first transition.')
catch ME
    help('models.triple_polytrope')
    rethrow(ME)
end

tof = TOFPlanet();

svec = linspace(1, 1/N, N);
[~, tind] = min(abs(svec-x(7)));
assert(tind > 2,...
    'First transition too close to surface; first polytrope has zero layers.')
[~, cind] = min(abs(svec-x(8)));
assert(cind > tind,...
    'Transitions are too close together; second polytrope has zero layers.')
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
