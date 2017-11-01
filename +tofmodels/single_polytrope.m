function tof = single_polytrope(N, x)
%SINGLE_POLYTROPE The simplest toy model planet.
%    SINGLE_POLYTROPE(N, x) returns an N-point TOFPlanet object with a
%    barotropes.Polytrope with constant x(1) and index x(2). The default level
%    spacing is one of equal radius increments between s/s0=1 and s/s0=1/N.

try
    narginchk(2,2)
    validateattributes(N, {'numeric'}, {'positive', 'integer'}, '', 'N', 1)
    validateattributes(x, {'numeric'}, {'vector', 'numel', 2, 'nonnegative'}, 2)
catch ME
    help('tofmodels.single_polytrope')
    rethrow(ME)
end

tof = TOFPlanet();
tof.si = linspace(1, 1/N, N);
tof.rhoi = ones(N,1);

eos1 = barotropes.Polytrope(x(1), x(2));
tof.eos = eos1;

end
