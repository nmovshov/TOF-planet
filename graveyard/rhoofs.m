function tof = rhoofs(svec, rhovec)
%RHOOFS Directly set level radii and densities for TOFPlanet.
%   RHOOFS(svec, rhovec) returns a TOFPlanet object whose si and rhoi properties
%   are initilized to input arguments svec and rhovec, respectively. The length of
%   svec and rhovec (which of course must be vectors of equal length) will
%   determine the number of levels in the model.
%
%   Obviously this can be easily done in any script or in the Command Window but
%   some driver scripts demand a handle to a parametrized model-building
%   function, so we supply one.
%
%   RHOOFS(srho) is a single input variant using an n-by-2 array instead of
%   two vectors. This simplifies calls to some function-functions, like fminsearch
%   and mhsample. The first column is svec and the second column is rhovec.

if nargin == 0, help('generators.rhoofs'), return, end
narginchk(1,2)
if nargin == 1
    validateattributes(svec,{'numeric'},{'2d','ncols',2},1)
    rhovec = svec(:,2);
    svec = svec(:,1);
else
    validateattributes(svec,{'numeric'},{'vector'},1)
    validateattributes(rhovec,{'numeric'},{'vector'},2)
    assert(length(svec) == length(rhovec), 'length(svec) ~= length(rhovec)')
end

tof = TOFPlanet();
tof.si = svec;
tof.rhoi = rhovec;

end
