function zvec = equal_volume(N)
%EQUAL_VOLUME Return a zvec distribution making layers of equal volume.
%    zvec = EQUAL_VOLUME(N) returns an N-vector of normalized radii such that
%    the volume of the spherical layer between zvec(i) and zvec(i+1) is constant.

validateattributes(N,{'numeric'},{'positive','integer','scalar'})

zvec = ((N:-1:1)/N).^(1/3);

end
