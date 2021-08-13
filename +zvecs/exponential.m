function zvec = exponential(N)
%EXPONENTIAL Return a zvec distribution with exponentially spaced radii.
%    zvec = EXPONENTIAL(N) returns an N-vector of exponentially spaced
%    normalized radii.

validateattributes(N,{'numeric'},{'positive','integer','scalar'})

zvec = log((exp(2)-1)*(N:-1:1)/N+1)/2;

end
