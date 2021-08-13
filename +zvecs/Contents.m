% +ZVECS TOF level distributions.
%
% Each function in this package returns a vector of length N suitable for use
% as the zvec argument to tof<n>.m or, more typically, to create the si
% property of a TOFPlanet object. For example:
%    tof = TOFPlanet();
%    tof.radius = 71492*1e3;
%    tof.si = tof.radius*zvecs.topheavy(1024);
%
% Distributions:
%   best         - Return our current notion of the best zvec distribution.
%   cosine       - Return a zvec distribution with cosine-like spacing.
%   equal_dr     - Return a zvec distribution with equally spaced radii.
%   equal_volume - Return a zvec distribution making layers of equal volume.
%   exponential  - Return a zvec distribution with exponentially spaced radii.
%   topheavy     - Return a zvec distribution with top-heavy spacing.
%   toptopheavy  - Return a zvec distribution with top-top-heavy spacing.
%   trizone      - Return a 3-zone zvec distribution.
