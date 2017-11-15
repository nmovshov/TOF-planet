function zvec = topheavy(N, skew)
%TOPHEAVY Return a zvec distribution with top-heavy spacing.
%    zvec = TOPHEAVY(N) returns a vector of normalized radii with three quarters
%    of them equally distributed in the top half and the rest equally distributed
%    in the bottom half of the interval (0,1].
%
%    zvec = TOPHEAVY(N, skew) returns a vector of normalized radii with skew(1) of
%    them equally distributed in the interval [1-skew(2), 1] and the rest equally
%    distributed in the interval (0, 1-skew(2)].

narginchk(1,2)
if (nargin < 2) || isempty(skew), skew = [3/4, 1/2]; end
validateattributes(N,{'numeric'},{'positive','integer','scalar'})
validateattributes(skew,{'numeric'},{'vector','numel',2,'>',0,'<',1})

n1 = fix(skew(1)*N);
n2 = N - n1;
lam1 = linspace(1, 1 - skew(2), n1);
lam2 = linspace(1 - skew(2), 1/N, n2+1);
zvec = flip(unique([lam1, lam2]));

end
