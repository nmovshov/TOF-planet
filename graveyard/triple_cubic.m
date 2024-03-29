function tof = triple_cubic(N, x, zstrat, forcemono)
%TRIPLE_CUBIC Piecewise-cubic density profile.
%    TRIPLE_CUBIC(N, x) returns an N-point TOFPlanet object with rhoi
%    approximated by a piecewise-cubic function:
%    
%        rho = x(1)*s^3 + x(2)*s^2  + x(3)*s  + x(4), if s > x(12),
%        rho = x(5)*s^3 + x(6)*s^2  + x(7)*s  + x(8), if x(13) < s < x(12),
%        rho = x(9)*s^3 + x(10)*s^2 + x(11)*s + 1,    otherwise.
%    
%    The default level spacing is one of equal radius increments between s/s0=1
%    and s/s0=1/N.
%    
%    TRIPLE_CUBIC(N, x, zstrat) lets you specify the zvec distribution
%    strategy. Pass a handle to a function that takes a single scalar integer
%    (number of layers) and returns a vector of that length with values in the
%    interval (0, 1], for normalized mean level radii.
%
%    TRIPLE_CUBIC(...,forcemono) where forcemono==true forces the resulting
%    density profile to be monotonically nonincreasing. Default forcemono is true.

if nargin == 0
    help('generators.triple_cubic')
    return
end
narginchk(2,4)
if ((nargin < 3) || isempty(zstrat)), zstrat = @zvecs.best; end
if ((nargin < 4) || isempty(forcemono)), forcemono = true; end
validateattributes(N, {'numeric'},{'positive','integer'},'','N',1)
validateattributes(x, {'numeric'},{'real','vector','numel',13},2)
validateattributes(zstrat, {'function_handle'},{},'','zstrat',3)
validateattributes(forcemono,{'logical'},{'scalar'},'','forcemono',4)
assert(x(12)>0 && x(12)<1, 'Transition (normalized) radius must be in (0,1).')
assert(x(13)>0 && x(13)<1, 'Second transition radius must be in (0,1).')
assert(x(13) <= x(12), 'Second transition must come before first transition.')

tof = TOFPlanet();

zvec = zstrat(N);
assert(isnumeric(zvec) && isvector(zvec) && (numel(zvec) == N),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')
assert(all(zvec > 0) && all(zvec <= 1),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')

a1 = x(1); b1 = x(2);  c1 = x(3);  d1 = x(4);
a2 = x(5); b2 = x(6);  c2 = x(7);  d2 = x(8);
a3 = x(9); b3 = x(10); c3 = x(11); d3 = 1;
zt = x(12);
zc = x(13);

dvec = NaN(1,N);
for k=1:N
    if zvec(k) > zt
        dvec(k) = a1*zvec(k)^3 + b1*zvec(k)^2 + c1*zvec(k) + d1;
    elseif zvec(k) > zc
        dvec(k) = a2*zvec(k)^3 + b2*zvec(k)^2 + c2*zvec(k) + d2;
    else
        dvec(k) = a3*zvec(k)^3 + b3*zvec(k)^2 + c3*zvec(k) + d3;
    end
end

if forcemono
    dvec(1) = max(dvec(1), 0);
    for k=2:N
        dvec(k) = max(dvec(k), dvec(k-1));
    end
end

tof.si = zvec;
tof.rhoi = dvec;

end
