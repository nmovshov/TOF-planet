function tof = tqrerep(N, x, zstrat, forcemono)
%TQREREP Piecewise-quadratic density profile, reparameterized...again.
%    TQREREP(N, x) returns an N-point TOFPlanet object with rhoi approximated by a
%    piecewise-quadratic function. Each quadratic segment is defined by its end
%    points and a curvature coefficient. The uppermost segment passes through the
%    point (1,0) and the lowermost passes through (0,1), leaving 9 free parameters
%    ordered as follows:
%    
%        x = [a1, y1, a2, ldy21, ldy22, a3, ldy32, rt, rc]
%    where
%        a1: curvature of first segment (a1*x^2 + b1*x + c1)
%        y1: rightlimit of rho(rt)
%        a2: curvature of second segment (a2*x^2 + b2*x + c2)
%        ldy21: log(leftlimit(rho(rt)) - rightlimit(rho(rt)))
%        ldy22: log(rightlimit(rho(rc)) - leftlimit(rho(rt)))
%        a3: curvature of third segment (a3*x^2 + b3*x + c3)
%        ld32: log(leftilimit(rho(rc)) - rightlimit(rho(rc)))
%        rt: location of first (upper) break point
%        rc: location of second (lower) break point
%    
%    The default level spacing is one of equal radius increments between s/s0=1
%    and s/s0=1/N.
%    
%    TQREREP(N, x, zstrat) lets you specify the zvec distribution strategy. Pass a
%    handle to a function that takes a single scalar integer (number of layers)
%    and returns a vector of that length with values in the interval (0, 1], for
%    normalized mean level radii.
%
%    TQREREP(...,forcemono) where forcemono==true forces the resulting density
%    profile to be monotonically nonincreasing. Default forcemono is false.

if nargin == 0 && nargout == 0
    help('tofmodels.tqrerep')
    return
end
narginchk(2,4)
if ((nargin < 3) || isempty(zstrat)), zstrat = @zvecs.best; end
if ((nargin < 4) || isempty(forcemono)), forcemono = false; end
validateattributes(N,{'numeric'},{'positive','integer'},'','N',1)
validateattributes(x,{'numeric'},{'real','vector','numel',9},2)
validateattributes(zstrat,{'function_handle'},{},'','zstrat',3)
validateattributes(forcemono,{'logical'},{'scalar'},'','forcemono',4)
assert(x(8)>0 && x(8)<1, 'Transition (normalized) radius must be in (0,1).')
assert(x(9)>0 && x(9)<1, 'Second transition radius must be in (0,1).')
assert(x(9) <= x(8), 'Second transition must come before first transition.')

tof = TOFPlanet();

zvec = zstrat(N);
assert(isnumeric(zvec) && isvector(zvec) && (numel(zvec) == N),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')
assert(all(zvec > 0) && all(zvec <= 1),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')

a1 = x(1); rot1 = x(2);
a2 = x(3); rot2 = rot1 + exp(x(4)); roc2 = rot2 + exp(x(5));
a3 = x(6); roc3 = roc2 + exp(x(7));
rt = x(8); rc = x(9);

% Upper envelope region
b1 = rot1/(rt - 1) - a1*(rt + 1);
c1 = -a1 - b1;

% Lower envelope region
b2 = (roc2 - rot2)/(rc - rt) - a2*(rc + rt);
c2 = rot2 - a2*rt^2 - b2*rt;

% Core region
b3 = (roc3 - 1)/rc - a3*rc;
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
