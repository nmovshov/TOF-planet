function tof = double_polytrope_w_core(N, x, zstrat)
%DOUBLE_POLYTROPE_W_CORE Model planet with two polytropes and a solid core.
%    DOUBLE_POLYTROPE_W_CORE(N, x) returns an N-point TOFPlanet object with
%    three barotropes eos objects. The first is a polytrope defined by constant
%    x(1) and index x(2) and applied to levels 1:tind-1. The second is a
%    polytrope defined by constant x(3) and index x(4) and applied to levels
%    tind:cind-1. The third is a constant density solid defined by x(5) and
%    applied to levels cind:N. Transition from first to second polytrope is at
%    level tind, the level with (si/s0) nearest x(6). Transition from second
%    polytrope to solid is at level cind, the level with (si/s0) nearest x(7).
%
%    The default level spacing is one of equal radius increments between s/s0=1
%    and s/s0=1/N, except that the transition levels are "snapped" to the given
%    values, x(6) and x(7).
%
%    DOUBLE_POLYTROPE_W_CORE(N, x, zstrat) lets you specify the level spacing.
%    Pass a handle to a function that takes a single scalar integer (number of
%    layers) and returns a vector of that length with values in the interval
%    (0,1], for normalized level radii. See +zvecs for examples.

if nargin == 0
    help('generators.double_polytrope_w_core')
    return
end
narginchk(2,3)
if ((nargin < 3) || isempty(zstrat)), zstrat = @(n)linspace(1, 1/n, n); end
validateattributes(N,{'numeric'},{'positive','integer'},'','N',1)
validateattributes(x,{'numeric'},{'nonnegative','vector','numel',7},2)
validateattributes(zstrat,{'function_handle'},{},'','zstrat',3)
assert(x(6) > 0 && x(6) < 1, 'Envelope transition radius must be in (0,1).')
assert(x(7) > 0 && x(7) < 1, 'Core radius must be in (0,1).')
assert(x(7) <= x(6), 'Core must be interior to lower envelope.')

% Generate and verify the layer distribution
zvec = zstrat(N);
assert(isnumeric(zvec) && isvector(zvec) && (numel(zvec) == N),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')
assert(all(zvec > 0) && all(zvec <= 1),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')

% Find and verify eos transitions
[~, tind] = min(abs(zvec-x(6)));
[~, cind] = min(abs(zvec-x(7)));
if tind == 1
    warning('TOF:GEN',...
        'Transition too close to surface; first polytrope will have zero levels.')
end
if tind == N
    warning('TOF:GEN',...
        'Transition too close to core; second polytrope will have zero levels.')
end

% Snap-to-grid transition locations
zvec(tind) = x(6);
zvec(cind) = x(7);

% Create the TOFPlanet object with assigned level distribution
tof = TOFPlanet();
tof.si = zvec;
tof.rhoi = ones(N,1);

% Create and assign the eoss
eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.Polytrope(x(3), x(4));
eos3 = barotropes.ConstDensity(x(5));

tof.eos = [repmat(eos1, tind - 1, 1);...
           repmat(eos2, cind - tind, 1);...
           repmat(eos3, N - cind + 1, 1)];

% Initialize density at 1-bar values, just because
for k=1:N
    tof.rhoi(k) = tof.eos(k).density(1e5);
end

end
