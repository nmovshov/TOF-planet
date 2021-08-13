function tof = layered_planet(N, EOSs, rts, zstrat, forcematch)
%LAYERED_PLANET Convenience generator for assigining EOSs to layers.
%    LAYERED_PLANET(N, EOSs, rts) returns an N-point TOFPlanet object with
%    length(EOSs) eos objects. Transition from EOSs(i) to EOSs(i+1) is at level
%    tind(i), the level with (si/s0) nearest rts(i). EOSs(1) is applied to
%    levels 1:tind(1)-1, EOSs(2) is applied to levels tind(1):tind(2)-1, and so
%    forth. EOSs(end) is applied to levels tind(end):N. Thus length(rts) must
%    equal length(EOSs)-1. The default level spacing is one of equal radius
%    increments between s/s0=1 and s/s0=1/N, except that the transition levels
%    are "snapped" to the given values in rts.
%
%    LAYERED_PLANET(N, EOSs, rts, zstrat) lets you specify the level-radii
%    distribution. Pass a handle to a function that takes a single scalar
%    integer (number of levels) and returns a vector of that length with values
%    in the interval (0, 1], for normalized mean level radii. Note that some
%    level radii may be slightly modified to match the requested transition
%    radii. A collection of pre-made distributions is available in package
%    +zvecs.
%
%    LAYERED_PLANET(...,forcematch) if forcematch=true forces the normalized
%    radii of the transition from EOSs(i) to EOSs(i+1) to exactly match rts(i).
%    This is applied _after_ the initial radii spacing. The default is
%    forcematch=true.
%
% Examples:
%    % Make a planet with a 2-polytrope envelope and a constant density core
%    eoss = [barotropes.Polytrope(1e5, 1),...
%            barotropes.Polytrope(8e4, 1),...
%            barotropes.ConstDensity(8000)];
%    tof = generators.layered_planet(128, eoss, [0.8, 0.2]);

% Handle inputs
if nargin == 0, help('generators.layered_planet'), return, end
narginchk(3,5)
if ((nargin < 4) || isempty(zstrat)), zstrat = @(n)linspace(1,1/n,n); end
if ((nargin < 5) || isempty(forcematch)), forcematch = true; end
validateattributes(N,{'numeric'},{'positive','integer'},'','N',1)
validateattributes(EOSs,{'barotropes.Barotrope'},{'vector'},'','EOSs',2)
validateattributes(rts,{'numeric'},{'vector'},'','rts',3)
validateattributes(zstrat,{'function_handle'},{},'','zstrat',4)
validateattributes(forcematch,{'logical'},{'scalar'},'','forcematch',5)
assert(all(rts > 0 & rts < 1), 'Transition radii must be in (0,1).')
assert(isequal(rts, sort(rts,'descend')), 'Levels are ordered from top down.')
assert(length(rts) == length(EOSs) - 1, 'Specify n-1 radii for n EOSs.')

% Generate and verify the level-radii distribution
zvec = zstrat(N);
assert(isnumeric(zvec) && isvector(zvec) && (numel(zvec) == N),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')
assert(all(zvec > 0) && all(zvec <= 1),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')

% Find and verify eos transition(s)
tind = nan(size(rts));
for k=1:length(rts)
    [~, tind(k)] = min(abs(zvec-rts(k)));
end
if tind(1) == 1
    warning('TOF:GEN',...
        'First transition too close to surface; eos 1 will have zero levels.')
end
for k=1:length(tind)-1
    if tind(k+1) == tind(k)
        warning('TOF:GEN',...
            'Transitions too close together; eos(%g) will have zero levels.', k)
    end
end
if tind(end) == N
    warning('TOF:GEN',...
        'Last transition too close to center; last eos will have just one level.')
end

% Snap-to-grid transition locations
if forcematch, zvec(tind) = rts; end

% Create the TOFPlanet object with assigned level-radii distribution
tof = TOFPlanet;
tof.si = zvec;
tof.rhoi = ones(N,1);

% Assign the eoss
tof.eos = repmat(barotropes.ConstDensity(0), N, 1);
for k=1:N
    zind = find(rts < zvec(k), 1);
    if ~isempty(zind)
        tof.eos(k) = EOSs(zind);
    else
        tof.eos(k) = EOSs(end);
    end
end

% Initialize a 1-bar density throughout (helps visualize the eos layers)
for k=1:N
    tof.rhoi(k) = tof.eos(k).density(1e5);
end

end
