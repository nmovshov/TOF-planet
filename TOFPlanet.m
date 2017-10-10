classdef TOFPlanet < handle
    %TOFPLANET Interior model of rotating fluid planet.
    %   This class implements a model of a rotating fluid planet using Theory of
    %   Figures to calculate the hydrostatic equilibrium shape and resulting
    %   gravity field. A TOFPlanet object is defined by a given mass, mean radius,
    %   rotation parameter, and optionally a barotrope.
    
    %% Properties
    properties
        name % model name
        Rm   % mean radius
        M    % mass
        mrot % rotation parameter, w^2s0^3/GM
        eos  % barotrope(s) (tip: help barotropes for options)
        opts % holds user configurable options
    end
    properties (Dependent)
        si   % vector of mean radii, top down
        rhoi % vector of densities on si grid
    end
    properties
        zvec % mean radii normalized by s0
        dvec % density normalized by mean density
    end
    properties (SetAccess = private)
        ss   % shape functions (returned by tof4.m)
        SS   % shape functions (returned by tof4.m)
        Js   % external gravity coefficients (returned by tof4.m)
    end
    properties (Dependent)
        M_calc
        rhobar
        qrot
        Req
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
            obj.zvec = linspace(0,1,N)';
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
                fprintf('  Relaxing to hydrostatic equilibrium...\n')
            end
            
            tic
            [obj.Js, out] = tof4(obj.zvec, obj.dvec, obj.mrot, obj.opts.dJtol);
            ET = toc;
            obj.ss = out.ss;
            obj.SS = out.SS;
            obj.aos = out.a0;
            
            if (obj.opts.verbosity > 1)
                fprintf('  Relaxing to hydrostatic equilibrium...done.\n')
                fprintf('  Elapsed time %s\n',lower(seconds2human(ET)))
            end
        end
        
        function r = level_surface(obj, l, theta)
            % Return r_l(theta) by expansion of Legendre polynomials.
            %
            % Usage:r = tof(l, theta)
            
            validateattributes(l,{'numeric'},{'positive','scalar','<=',1},'','l',1)
            validateattributes(theta,{'numeric'},{'vector','>=',0,'<=',2*pi},'','theta',2)
            
            k = find(obj.zvec >= l, 1, 'first');
            mu = cos(theta);
            shp = obj.ss.s0(k)*Pn(0,mu) + obj.ss.s2(k)*Pn(2,mu) + ...
                obj.ss.s4(k)*Pn(4,mu) + obj.ss.s6(k)*Pn(6,mu) + ...
                obj.ss.s8(k)*Pn(8,mu);
            r = obj.zvec(k)*(1 + shp);
        end
    end % End of public methods block
    
    %% Private methods
    methods (Access = private)
    end % End of private methods block
    
    %% Access methods
    methods
        function val = get.qrot(obj)
            val = obj.mrot*obj.aos^3;
        end
        
        function val = get.Req(obj)
            val = obj.aos*obj.Rm;
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

%% Class-related functions
function y = Pn(n,x)
% Fast implementation of ordinary Legendre polynomials of low degree.
switch n
    case 0
        y = ones(size(x));
    case 1
        y = x;
    case 2
        y = 0.5*(3*x.^2 - 1);
    case 3
        y = 0.5*(5*x.^3 - 3*x);
    case 4
        y = (1/8)*(35*x.^4 - 30*x.^2 + 3);
    case 5
        y = (1/8)*(63*x.^5 - 70*x.^3 + 15*x);
    case 6
        y = (1/16)*(231*x.^6 - 315*x.^4 + 105*x.^2 - 5);
    case 7
        y = (1/16)*(429*x.^7 - 693*x.^5 + 315*x.^3 - 35*x);
    case 8
        y = (1/128)*(6435*x.^8 - 12012*x.^6 + 6930*x.^4 - 1260*x.^2 + 35);
    case 9
        y = (1/128)*(12155*x.^9 - 25740*x.^7 + 18018*x.^5 - 4620*x.^3 + 315*x);
    case 10
        y = (1/256)*(46189*x.^10 - 109395*x.^8 + 90090*x.^6 - 30030*x.^4 + 3465*x.^2 - 63);
    case 11
        y = (1/256)*(88179*x.^11 - 230945*x.^9 + 218790*x.^7 - 90090*x.^5 + 15015*x.^3 - 693*x);
    case 12
        y = (1/1024)*(676039*x.^12 - 1939938*x.^10 + 2078505*x.^8 - 1021020*x.^6 + 225225*x.^4 - 18018*x.^2 + 231);
    otherwise
        assert(isvector(x))
        Pnm = legendre(n,x);
        y = Pnm(1,:);
        if ~isrow(x), y = y'; end
end
end
