function tof = triple_polytrope(N, x, zstrat)
%TRIPLE_POLYTROPE Model planet approximated by three polytropes.
%    TRIPLE_POLYTROPE(N, x) returns an N-point TOFPlanet object with three
%    barotropes.Polytrope eos objects. First polytrope defined by constant x(1)
%    and index x(2) is applied to levels 1:tind-1. Second polytrope defined by
%    constant x(3) and index x(4) is applied to levels tind:cind-1. Third
%    polytrope defined by constant x(5) and index x(6) is applied to levels
%    cind:N. Transition from first to second polytrope is at level tind, the
%    level with (si/s0) nearest x(7). Transition from second to third polytrope
%    is at level cind, the level with (si/s0) nearest x(8). The default level
%    spacing is one of equal radius increments between s/s0=1 and s/s0=1/N,
%    except that the transition levels are snapped to the given values, x(7)
%    and x(8).
%
%    TRIPLE_POLYTROPE(N, x, zstrat) lets you specify the level-radii
%    distribution. Pass a handle to a function that takes a single scalar
%    integer (number of levels) and returns a vector of that length with values
%    in the interval (0, 1], for normalized mean level radii. For example, to
%    set levels with equally spaced radii use zstrat=@(n)linspace(1,1/n,n).
%    Note that some level radii may be slightly modified to match the requested
%    transition radii. A collection of pre-made distributions is available in
%    package +zvecs.

if nargin == 0, help('generators.triple_polytrope'), return, end
narginchk(2,3)
if ((nargin < 3) || isempty(zstrat)), zstrat = @(n)linspace(1,1/n,n); end
validateattributes(N,{'numeric'},{'positive','integer'},'','N',1)
validateattributes(x,{'numeric'},{'vector','numel',8,'nonnegative'},2)
validateattributes(zstrat,{'function_handle'},{},'','zstrat',3)
assert(x(7) > 0 && x(7) < 1, 'First transition radius must be in (0,1).')
assert(x(8) > 0 && x(8) < 1, 'Second transition radius must be in (0,1).')
assert(x(8) <= x(7), 'Second transition should be deeper than the first.')

% Generate and verify the layer distribution
zvec = zstrat(N);
assert(isnumeric(zvec) && isvector(zvec) && (numel(zvec) == N),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')
assert(all(zvec > 0) && all(zvec <= 1),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')

% Find and verify eos transitions
[~, tind] = min(abs(zvec-x(7)));
[~, cind] = min(abs(zvec-x(8)));
if tind == 1
    warning('TOF:GEN',...
        'First transition too close to surface; first polytrope has zero levels.')
end
if cind == tind
    warning('TOF:GEN',...
        'Transitions are too close together; second polytrope has zero levels.')
end
if cind == N
    warning('TOF:GEN',...
        'Second transition too deep; last polytrope has just one level.')
end

% Snap-to-grid transition locations
zvec(tind) = x(7);
zvec(cind) = x(8);

% Create the TOFPlanet object with assigned layer distribution
tof = TOFPlanet;
tof.si = zvec;
tof.rhoi = ones(N,1);

% Create and assign the eoss
eos1 = barotropes.Polytrope(x(1), x(2));
eos2 = barotropes.Polytrope(x(3), x(4));
eos3 = barotropes.Polytrope(x(5), x(6));
tof.eos = [repmat(eos1, tind - 1, 1);...
           repmat(eos2, cind - tind, 1);...
           repmat(eos3, N - cind + 1, 1)];

% Initialize density at 1-bar values, just because
for k=1:N
    tof.rhoi(k) = tof.eos(k).density(1e5);
end

end
