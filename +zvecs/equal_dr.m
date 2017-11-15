function zvec = equal_dr(N)
%EQUAL_DR Return a zvec distribution with equally spaced radii.
%    zvec = EQUAL_DR(N) returns an N-vector of equally spaced normalized radii. So
%    really it's another name for linapce(1,1/N,N) but we include it for
%    completeness.

validateattributes(N,{'numeric'},{'positive','integer','scalar'})

zvec = linspace(1, 1/N, N);

end
