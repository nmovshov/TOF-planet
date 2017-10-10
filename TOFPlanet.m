classdef TOFPlanet < handle
    %TOFPLANET Interior model of rotating fluid planet.
    %   This class implements a model of a rotating fluid planet using Theory of
    %   Figures to calculate the hydrostatic equilibrium shape and resulting
    %   gravity field. A TOFPlanet object is defined by a given mass, mean radius,
    %   rotation parameter, and optionally a barotrope.
    
    %% Properties
    properties
        name % model name
        s0   % mean radius
        M    % mass
        mrot % rotation parameter, w^2s0^3/GM
        eos  % barotrope(s) (tip: help barotropes for options)
        
        zvec % mean radii normalized by s0
        dvec % density normalized by mean density
        
        opts % holds user configurable options
    end
    
    %% A simple constructor
    methods
        function obj = TOFPlanet(N, varargin)
            % A simple constructor of TOFPlanet objects.
            % TOFPlanet(N, 'OPTION1', VALUE, 'OPTION2', VALUE2,...)
            
            % The grid size os required
            if nargin == 0
                error(['Required argument missing: specify number of radius',...
                    ' grid points as first input argument.'])
            end
            validateattributes(N,{'numeric'},{'positive','integer','scalar'},...
                '','N',1)
            
            % Populate options struct
            obj.opts = tofset(varargin{:});
            
            % Set a default density grid - usually replaced by user
            obj.zvec = linspace(1,0,N)';
            obj.dvec = ones(size(obj.zvec));
            
        end
    end % End of constructor block
    
    %% Public methods
    methods (Access = public)
    end % End of public methods block
    
    %% Private methods
    methods (Access = private)
    end % End of private methods block
    
    %% Access methods
    methods
    end % End of access methods block
    
    %% Static methods
    methods (Static)
    end % End of static methods block
end % End of classdef

