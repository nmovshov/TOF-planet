function [zvec,dvec] = ppwd_profile(zvec,x,rho0,forcemono)
%PPWD_PROFILE Polynomial planet with smoothed-step discontinuities.
%    prof = PPWD_PROFILE(N, x, rho0) returns an N-level density profile
%    approximated by a single polynomial of normalized mean radius z=s/s0
%    overlain with two smoothed-step-functions. The default radii spacing is
%    one of equal increments between z=1 and z=1/N. The parameter vector x and
%    the surface density rho0 together fully define the density profile.
%
%    The first three elements of x define the location, scale, and sharpness of
%    the first (inner) step. Specify location 0<z_1<1 in normalized mean radius
%    and specify scale in kg/m^3. The sharpness parameter is non-dimensional
%    positive; a larger value results in a narrower step. Try ~100. The next
%    three elements of x specify the same parameters for the second (outer)
%    step. The remaining elements are the coefficients of the polynomial.
%    Remember that a deg-n polynomial will be constrained by boundary
%    conditions and defined with n-1 coefficients rather than n+1.
%
%    PPWD_PROFILE(zvec,...) specifies the level radii distribution directly,
%    rather than a number of levels. The vector zvec should have normalized
%    values in (0,1] and sorted top-down (rho0 is applied and zvec(1)).
%
%    PPWD_PROFILE(...,forcemono) where forcemono==true forces the resulting
%    density profile to be monotonically non increasing. This shouldn't be
%    necessary but occasionally for high-resolution models and high-order
%    polynomials there is an unwanted min somewhere in the first few layers.

if nargin == 0
    fprintf('Usage:\n\tprof=ppwd_profile(N, x, rho0, [forcemono=true])\n')
    fprintf(' or\n')
    fprintf('\tprof=ppwd_profile(zvec, x, rho0, [forcemono=true])\n')
    fprintf('\tx(1): inner step location; 0 < z1=s/s0 < 1\n')
    fprintf('\tx(2): inner step height (in real density units)\n')
    fprintf('\tx(3): inner step dimensionless sharpness (try ~100)\n')
    fprintf('\tx(4): outer step location; 0 < z2=s/s0 < 1\n')
    fprintf('\tx(5): outer step height (in real density units)\n')
    fprintf('\tx(6): outer step dimensionless sharpness (try ~100)\n')
    fprintf('\tx(7:end): polynomial coefficients in decreasing order\n')
    return
end
narginchk(3,4)
if nargin < 4 || isempty(forcemono), forcemono = true; end
if isscalar(zvec) % hack to use N equally spaced levels
    N = zvec;
    zvec = linspace(1, 1/N, N);
end
assert(isnumeric(zvec) && isvector(zvec),...
    'zvec must be a vector with values in (0,1].')
assert(all(zvec > 0) && all(zvec <= 1),...
    'zvec must be a vector of length N with values in (0,1].')

% First, the polynomial
y = [x(7:end), 0, rho0 - polyval([x(7:end),0,0],zvec(1))];
dvec = polyval(y,zvec);
if forcemono
    dvec(1) = max(dvec(1), 0);
    for k=2:length(dvec)
        dvec(k) = max(dvec(k), dvec(k-1));
    end
end
dvec = dvec(:); % make sure it's a column
zvec = zvec(:); % make sure it's a column

% Add density jumps
ro0 = dvec(1);

z = x(4); scale = x(5); sharpness = x(6); % outer
y = scale*(pi/2 + atan(-sharpness*(zvec - z)))/pi;
dvec = dvec + y;

z = x(1); scale = x(2); sharpness = x(3); % inner
y = scale*(pi/2 + atan(-sharpness*(zvec - z)))/pi;
dvec = dvec + y;

dvec = dvec - dvec(1) + ro0;

% Undocumented convenience hack: return as 2 vectors if nargout>1
if nargout < 2, zvec = [zvec, dvec]; end
end
