function tof = tripquad_v1(N, x, forcemono, zstrat)
%TRIPQUAD_V1 Piecewise-quadratic density profile.
%    TRIPQUAD_V1(N, x) returns an N-point TOFPlanet object with normalized rhoi
%    (which we call D) approximated by a piecewise-quadratic function in
%    normalized radius (which we call Z). Each quadratic segment is defined by its
%    three polynomial coefficients (a_i, b_i, c_i). The segment breakpoints are at
%    z1 (upper) and z2 (lower). The uppermost segment passes through the point
%    (z,d)=(1,0) and the lowermost passes through (z,d)=(0,1) leaving 9 free
%    parameters. They are ordered as follows:
%        x = [a1, b1, a2, b2, c2, a3, b3, z1, z2]
%    where
%        d_i(z) = a_i*z^2 + b_i*z + c_i
%    and the boundary conditions require c1=-a1-b1 and c3=1.
%    
%    The default level spacing is one of equal radius increments between z=1 and
%    z=1/N.
%    
%    TRIPQUAD_V1(...,forcemono) where forcemono==true forces the resulting density
%    profile to be monotonically nonincreasing. Default forcemono is false.

if nargin == 0 && nargout == 0
    help('tofmodels.tripquad_v1')
    return
end
narginchk(2,4)
if ((nargin < 3) || isempty(forcemono)), forcemono = false; end
if ((nargin < 4) || isempty(zstrat)), zstrat = @zvecs.best; end
validateattributes(N,{'numeric'},{'positive','integer'},'','N',1)
validateattributes(x,{'numeric'},{'real','vector','numel',9},2)
validateattributes(forcemono,{'logical'},{'scalar'},'','forcemono',3)
validateattributes(zstrat,{'function_handle'},{},'','zstrat',4)
assert(x(8)>0 && x(8)<1, 'z1 must be in (0,1).')
assert(x(9)>0 && x(9)<1, 'z2 must be in (0,1).')
assert(x(9) < x(8), 'z2 must less than z1.')

tof = TOFPlanet();

zvec = zstrat(N);
assert(isnumeric(zvec) && isvector(zvec) && (numel(zvec) == N),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')
assert(all(zvec > 0) && all(zvec <= 1),...
    '@zstrat(N) must return a vector of length N with values in (0,1].')

a1 = x(1); b1 = x(2); c1 = -a1 - b1;
a2 = x(3); b2 = x(4); c2 = x(5);
a3 = x(6); b3 = x(7); c3 = 1;
z1 = x(8);
z2 = x(9);

dvec = NaN(1,N);
for k=1:N
    if zvec(k) > z1
        dvec(k) = a1*zvec(k)^2 + b1*zvec(k) + c1;
    elseif zvec(k) > z2
        dvec(k) = a2*zvec(k)^2 + b2*zvec(k) + c2;
    else
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
