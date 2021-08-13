function tof = polynomial(N, x, zstrat, forcemono)
%POLYNOMIAL A single polynomial in normalized mean radius.
%    POLYNOMIAL(N, x) returns an N-point TOFPlanet object with rhoi approximated
%    by a single polynomial of z=s/s0, with coefficients x. The default level
%    spacing is one of equal radius increments between s/s0=1 and s/s0=1/N.
%    
%    POLYNOMIAL(N, x, zstrat) lets you specify the level-radii distribution.
%    Pass a handle to a function that takes a single scalar integer (number of
%    levels) and returns a vector of that length with values in the interval
%    (0, 1], for normalized mean level radii.
%
%    POLYNOMIAL(...,forcemono) where forcemono==true forces the resulting
%    density profile to be monotonically nonincreasing. Default forcemono is
%    false, because if we mess up the polynomial coefficients it's best to find
%    out early.

if nargin == 0
    fprintf('Usage:\n\ttof = polynomial(N, x, zstrat, forcemono)\n')
    return
end
narginchk(2,4)
if nargin < 3 || isempty(zstrat), zstrat = @(n)linspace(1, 1/n, n); end
if nargin < 4 || isempty(forcemono), forcemono = false; end
validateattributes(N,{'numeric'},{'positive','integer'},'','N',1)
validateattributes(x,{'numeric'},{'real','finite','vector'},'','x',2)
validateattributes(zstrat,{'function_handle'},{},'','zstrat',3)
validateattributes(forcemono,{'logical'},{'scalar'},'','forcemono',4)

% Generate and verify the level-radii distribution
zvec = zstrat(N);
assert(isnumeric(zvec) && isvector(zvec) && (numel(zvec) == N),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')
assert(all(zvec > 0) && all(zvec <= 1),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')

% Generate the density profile
dvec = polyval(x,zvec);
if forcemono
    dvec(1) = max(dvec(1), 0);
    for k=2:N
        dvec(k) = max(dvec(k), dvec(k-1));
    end
end

% Create the TOFPlanet object and assign density profile
tof = TOFPlanet();
tof.si = zvec;
tof.rhoi = dvec;

end
