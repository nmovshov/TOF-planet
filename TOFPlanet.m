classdef TOFPlanet < handle
    %TOFPLANET Interior model of rotating fluid planet.
    %   This class implements a model of a rotating fluid planet using Theory of
    %   Figures to calculate the hydrostatic equilibrium shape and resulting
    %   gravity field. A TOFPlanet object is defined by a given mass, mean radius,
    %   rotation parameter, and optionally a barotrope.
    
    %% Properties
    properties
        name % model name
        opts % holds user configurable options
    end
    
    %% A simple constructor
    methods
        function obj = TOFPlanet(varargin)
            % The constructor only populates the options struct.
            obj.opts = tofset(varargin{:});
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

