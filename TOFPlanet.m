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
    properties (SetAccess = private)
        ss   % shape functions (returned by tof4.m)
        SS   % shape functions (returned by tof4.m)
        Js   % external gravity coefficients (returned by tof4.m)
    end
    properties (Dependent)
        qrot
        a0
        J2
        J4
        J6
        J8
    end
    properties (Access = private)
        aos = 1
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
            
            % Initialize density grid (usually replaced by user!)
            obj.zvec = linspace(1,0,N)';
            obj.dvec = ones(size(obj.zvec));
            
            % Initialize shape functions
            ss.s0(N,1)=0; ss.s2(N,1)=0; ss.s4(N,1)=0; ss.s6(N,1)=0; ss.s8(N,1)=0;
            SS.S0(N,1)=0; SS.S2(N,1)=0; SS.S4(N,1)=0; SS.S6(N,1)=0; SS.S8(N,1)=0;
            obj.ss = ss;
            obj.SS = SS;
            
            % Set other defaults
            obj.mrot = 0;
            
        end
    end % End of constructor block
    
    %% Public methods
    methods (Access = public)
        function ET = relax_to_HE(obj)
            % Call tof4 to obtain equilibrium shape and gravity.
            
            if (obj.opts.verbosity > 1)
                fprintf('  Calculating hydrostatic equilibrium...\n')
            end
            
            tic
            [obj.Js, out] = tof4(obj.zvec, obj.dvec, obj.mrot);
            ET = toc;
            obj.ss = out.ss;
            obj.SS = out.SS;
            obj.aos = out.a0;
            
            if (obj.opts.verbosity > 1)
                fprintf('  Calculating hydrostatic equilibrium...done.\n')
                fprintf('  Elapsed time %s\n',lower(seconds2human(ET)))
            end
        end
    end % End of public methods block
    
    %% Private methods
    methods (Access = private)
    end % End of private methods block
    
    %% Access methods
    methods
        function val = get.a0(obj)
            val = obj.aos*obj.s0;
        end
        
        function val = get.J2(obj)
            val = obj.Js(2);
        end
        
        function val = get.J4(obj)
            val = obj.Js(3);
        end
        
        function val = get.J6(obj)
            val = obj.Js(4);
        end
        
        function val = get.J8(obj)
            val = obj.Js(5);
        end
        
    end % End of access methods block
    
    %% Static methods
    methods (Static)
    end % End of static methods block
end % End of classdef

