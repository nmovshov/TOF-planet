function tof = tqrep(N, x, zstrat, forcemono)
%TQREP Piecewise-quadratic density profile, reparameterized.
%    TQREP(N, x) returns an N-point TOFPlanet object with normalized rhoi (which
%    we call D) approximated by a piecewise-quadratic function in normalized
%    radius (which we call Z). Each quadratic segment is defined by its end points
%    and a curvature coefficient. The segment breakpoints are at z1 (upper) and z2
%    (lower). The density at z1 is d11 (the right-limit) and d21 (the
%    left-limit). The density at z2 is d22 (the right-limit) and d32 (the
%    left-limit). The uppermost segment passes through the point (z,d)=(1,0) and
%    the lowermost passes through (z,d)=(0,1) leaving 9 free parameters. They are
%    ordered as follows:
%    
%        x = [a1, y11, a2, y21, y22, a3, y32, lz1, lz2]
%    where
%        a1: curvature of first segment (a1*x^2 + b1*x + c1)
%        y11: log(d11) - log(1 - d11)
%        a2: curvature of second segment (a2*x^2 + b2*x + c2)
%        y21: log(1 - d21) - log(d21 - d11)
%        y22: log(1 - d22) - log(d22 - d21)
%        a3: curvature of third segment (a3*x^2 + b3*x + c3)
%        y32: log(1 - d32) - log(d32 - d22)
%        lz1: log(z1) - log(1 - z1)
%        lz2: log(z2) - log(z1 - z2)
%    
%    The default level spacing is one of equal radius increments between z=1 and
%    z=1/N.
%    
%    TQREP(N, x, zstrat) lets you specify the zvec distribution strategy. Pass a
%    handle to a function that takes a single scalar integer (number of layers)
%    and returns a vector of that length with values in the interval (0, 1], for
%    normalized mean level radii.
%
%    TQREP(...,forcemono) where forcemono==true forces the resulting density
%    profile to be monotonically nonincreasing. Default forcemono is false.

if nargin == 0 && nargout == 0
    help('tofmodels.tqrep')
    return
end
narginchk(2,4)
if ((nargin < 3) || isempty(zstrat)), zstrat = @zvecs.best; end
if ((nargin < 4) || isempty(forcemono)), forcemono = false; end
validateattributes(N,{'numeric'},{'positive','integer'},'','N',1)
validateattributes(x,{'numeric'},{'real','vector','numel',9},2)
validateattributes(zstrat,{'function_handle'},{},'','zstrat',3)
validateattributes(forcemono,{'logical'},{'scalar'},'','forcemono',4)

tof = TOFPlanet();

zvec = zstrat(N);
assert(isnumeric(zvec) && isvector(zvec) && (numel(zvec) == N),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')
assert(all(zvec > 0) && all(zvec <= 1),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')

a1 = x(1); y11 = x(2);
a2 = x(3); y21 = x(4); y22 = x(5);
a3 = x(6); y32 = x(7);
lz1 = x(8); lz2 = x(9);

bexp = @(x)min(exp(x), realmax/10); % bounded exp to avoid overflow
d11 = bexp(y11)/(1 + bexp(y11));
d21 = (1 + d11*bexp(y21))/(1 + bexp(y21));
d22 = (1 + d21*bexp(y22))/(1 + bexp(y22));
d32 = (1 + d22*bexp(y32))/(1 + bexp(y32));
rt = bexp(lz1)/(1 + bexp(lz1));
rc = rt*bexp(lz2)/(1 + bexp(lz2));

% Upper envelope region
b1 = d11/(rt - 1) - a1*(rt + 1);
c1 = -a1 - b1;

% Lower envelope region
b2 = (d22 - d21)/(rc - rt) - a2*(rc + rt);
c2 = d21 - a2*rt^2 - b2*rt;

% Core region
b3 = (d32 - 1)/rc - a3*rc;
c3 = 1;

dvec = NaN(1,N);
for k=1:N
    if zvec(k) > rt
        % Upper envelope
        dvec(k) = a1*zvec(k)^2 + b1*zvec(k) + c1;
    elseif zvec(k) > rc
        % Lower envelope
        dvec(k) = a2*zvec(k)^2 + b2*zvec(k) + c2;
    else
        % Core
        dvec(k) = a3*zvec(k)^2 + b3*zvec(k) + c3;
    end
end

if forcemono
    for k=2:N
        dvec(k) = max(dvec(k), dvec(k-1));
    end
end

tof.si = zvec;
tof.rhoi = dvec;

end
