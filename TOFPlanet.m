classdef TOFPlanet < handle
    %TOFPLANET Interior model of rotating fluid planet.
    %   This class implements a model of a rotating fluid planet using Theory of
    %   Figures to calculate the hydrostatic equilibrium shape and resulting
    %   gravity field. A TOFPlanet object is defined by a given mass, mean radius,
    %   rotation parameter, and optionally a barotrope.
    
    %% Properties
    properties
        name   % model name
        mass   % reference mass
        radius % reference radius (equatorial!)
        P0     % reference pressure
        si     % vector of mean radii, top down
        rhoi   % vector of densities on si grid
        mrot   % rotation parameter, w^2s0^3/GM
        eos    % barotrope(s) (tip: help barotropes for options)
        opts   % holds user configurable options (tip: help tofset)
    end
    properties (SetAccess = private)
        ss     % shape functions (returned by tof4.m)
        SS     % shape functions (returned by tof4.m)
        Js     % external gravity coefficients (returned by tof4.m)
        betam  % mass renormalization factor (obj.mass/obj.M)
        alfar  % radius renormalization factor (obj.radius/obj.a0)
    end
    properties (Dependent)
        M      % calculated mass
        s0     % calculated mean radius (another name for obj.si(1))
        a0     % calculated equatorial radius
        rhobar % calculated mean density
        qrot   % rotation parameter referenced to a0
        Ui     % gravitational potential on grid
        Pi     % hydrostatic pressure on grid
        J2     % convenience alias to obj.Js(2)
        J4     % convenience alias to obj.Js(3)
        J6     % convenience alias to obj.Js(4)
        J8     % convenience alias to obj.Js(5)
    end
    properties (GetAccess = private)
        aos;   % calculated equatorial to mean radius ratio (from tof4.m)
        N;     % convenience name for length(obj.si)
        G;     % Gravitational constant
        u;     % let's hold a units struct for convenience
    end
    
    %% A simple constructor
    methods
        function obj = TOFPlanet(N, varargin)
            % A simple constructor of TOFPlanet objects.
            % TOFPlanet(N, 'OPTION1', VALUE, 'OPTION2', VALUE2,...)
            
            % The grid size is required
            if nargin == 0
                error(['Required argument missing: specify number of radius',...
                    ' grid points as first input argument.'])
            end
            validateattributes(N,{'numeric'},{'positive','integer','scalar'},...
                '','N',1)
            
            % Populate options struct
            obj.opts = tofset(varargin{:});
            
            % Init privates
            obj.N = N;
            obj.aos = [];
            if obj.opts.debug
                obj.u = setUnits;
            else
                obj.u = setFUnits;
            end
            obj.G = obj.u.gravity;
        end
    end % End of constructor block
    
    %% Public methods
    methods (Access = public)
        function ET = relax_to_barotrope(obj)
            % Relax equilibrium shape functions, Js, and density simultaneously.
            
            if isempty(obj.eos)
                warning('Set valid barotrope first (<obj>.eos = <barotrope>)')
                return
            end
            if isempty(obj.P0)
                warning('Setting reference pressure to zero (<obj>.P0=0).')
                obj.P0 = 0*obj.u.bar;
            end
            
            t_rlx = tic;
            
            % Optional communication
            verb = obj.opts.verbosity;
            if (verb > 0)
                fprintf('Relaxing to desired barotrope...\n\n')
            end
            
            % Main loop
            iter = 1;
            while (iter <= obj.opts.MaxIterBar)
                t_pass = tic;
                if (verb > 0)
                    fprintf('Baropass %d (of max %d)...\n',...
                        iter, obj.opts.MaxIterBar)
                end
                
                old_Js = obj.Js;
                if isempty(old_Js), old_Js = [-1,0,0,0,0]; end
                old_ro = obj.rhoi;
                obj.relax_to_HE();
                obj.update_densities;
                dJs = abs((obj.Js - old_Js)./old_Js);
                dJs = max(double(dJs(isfinite(dJs))));
                dro = obj.rhoi./old_ro;
                dro = var(double(dro(isfinite(dro))));
                
                if (verb > 0)
                    fprintf('Baropass %d (of max %d)...done. (%g sec)\n',...
                        iter, obj.opts.MaxIterBar, toc(t_pass))
                    fprintf('var(drho) = %g (%g required); dJ = %g (%g required).\n\n',...
                        dro, obj.opts.drhotol, dJs, obj.opts.dJtol)
                end
                
                % The stopping criterion is to satisfy both J and rho tolerance
                if (dro < obj.opts.drhotol) && all(dJs < obj.opts.dJtol)
                    break
                end
                
                % end the main loop
                iter = iter + 1;
            end
            ET = toc(t_rlx);
            
            % Renormalize densities, radii, Js.
            obj.alfar = obj.radius/obj.a0;
            obj.si = obj.si*obj.alfar;
            obj.betam = obj.mass/obj.M;
            obj.rhoi = obj.rhoi*obj.betam;
            
            % Optional communication
            if (verb > 0)
                fprintf('Relaxing to desired barotrope...done.\n')
                fprintf('Total elapsed time %s\n',lower(seconds2human(ET)))
            end
        end
        
        function [ET, dJ] = relax_to_HE(obj)
            % Call tof4 to obtain equilibrium shape and gravity.
            
            if (obj.opts.verbosity > 1)
                fprintf('  Relaxing to hydrostatic equilibrium...\n')
            end
            
            tic
            zvec = double(obj.si/obj.si(1));
            dvec = double(obj.rhoi/obj.rhobar);
            [obj.Js, out] = tof4(zvec, dvec, obj.mrot, obj.opts.dJtol, obj.opts.MaxIterHE);
            ET = toc;
            dJ = out.dJs;
            obj.ss = structfun(@flipud, out.ss, 'UniformOutput', false);
            obj.SS = structfun(@flipud, out.SS, 'UniformOutput', false);
            obj.aos = out.a0;
            
            if (obj.opts.verbosity > 1)
                fprintf('  Relaxing to hydrostatic equilibrium...done.\n')
                fprintf('  Elapsed time %s\n',lower(seconds2human(ET)))
            end
        end
        
        function dro = update_densities(obj)
            % Set level surface densities to match prescribed barotrope.
            
            if isempty(obj.eos)
                warning('Make sure input barotrope (<obj>.eos) is set.')
                return
            end
            
            t_rho = tic;
            verb = obj.opts.verbosity;
            if (verb > 1)
                fprintf('  Updating level surface densities...')
            end
            P = obj.Pi;
            if isscalar(obj.eos)
                newro = obj.eos.density(P);
            else
                newro = repmat(obj.rhoi(1), obj.N, 1);
                for k=1:length(newro)
                    newro(k) = obj.eos(k).density(P(k));
                end
            end
            dro = ((newro - obj.rhoi)./obj.rhoi);
            if (verb > 2)
                fprintf('done. (%g sec)\n', toc(t_rho))
            elseif (verb > 1)
                fprintf('done.\n')
            end
            obj.rhoi = newro;
        end
        
        function r = level_surface(obj, el, theta)
            % Return r_l(theta) by expansion of Legendre polynomials.
            %
            % Usage:r = tof(el, theta)
            
            validateattributes(el,{'numeric'},{'positive','scalar','<=',1},'','l',1)
            validateattributes(theta,{'numeric'},{'vector','>=',0,'<=',2*pi},'','theta',2)
            
            k = find(obj.si/obj.si(1) <= el, 1, 'first');
            mu = cos(theta);
            shp = obj.ss.s0(k)*Pn(0,mu) + obj.ss.s2(k)*Pn(2,mu) + ...
                obj.ss.s4(k)*Pn(4,mu) + obj.ss.s6(k)*Pn(6,mu) + ...
                obj.ss.s8(k)*Pn(8,mu);
            r = obj.si(k)*(1 + shp);
        end
        
        function m = total_mass(obj,method)
            if nargin == 1, method = 'layerz'; end
            switch lower(method)
                case 'trapz'
                    m = -4*pi*trapz(double(obj.si), double(obj.rhoi.*obj.si.^2));
                    m = m*obj.u.kg;
                case 'layerz'
                    drho = [obj.rhoi(1); diff(obj.rhoi)];
                    m = (4*pi/3)*sum(drho.*obj.si.^3);
                case 'integralz'
                    x = double(obj.si);
                    v = double(obj.rhoi.*obj.si.^2);
                    fun = @(z)interp1(x, v, z, 'pchip');
                    m = 4*pi*integral(fun, 0 , x(1));
                otherwise
                    error('Unknown mass calculation method.')
            end
        end
        
        function [ah, lh, gh] = plot_barotrope(obj, varargin)
            % Plot P(rho) of current model and optionally of input barotrope.
            
            % Don't bother if there is no pressure
            if isempty(obj.Pi)
                warning('Uninitialized object.')
                return
            end
            
            % Input parsing
            p = inputParser;
            p.FunctionName = mfilename;
            p.addParameter('axes', [],...
                @(x)isscalar(x) && isgraphics(x,'axes') && isvalid(x));
            p.addParameter('showinput', false,...
                @(x)isscalar(x) && islogical(x));
            p.addParameter('showscaledinput', false,...
                @(x)isscalar(x) && islogical(x));
            p.addParameter('includecore', false,...
                @(x)isscalar(x) && islogical(x));
            p.parse(varargin{:})
            pr = p.Results;
            
            % Prepare the canvas
            if isempty(pr.axes)
                fh = figure;
                set(fh, 'defaultTextInterpreter', 'latex')
                set(fh, 'defaultLegendInterpreter', 'latex')
                ah = axes;
                hold(ah, 'on')
            else
                ah = pr.axes;
                hold(ah, 'on')
            end
            
            % Prepare the data: model
            x_tof = double(obj.rhoi);
            y_tof = double(obj.Pi);
            
            % Prepare the data: input
            if pr.showinput && ~isempty(obj.eos) && (range(x_tof) > 0)
                x_bar = linspace(min(x_tof), max(x_tof));
                if isscalar(obj.eos)
                    y_bar = double(obj.eos.pressure(x_bar));
                else
                    v = 1:length(unique(x_tof));
                    ind = interp1(unique(x_tof), v, x_bar, 'nearest', 'extrap');
                    y_bar = nan(size(x_bar));
                    for k=1:length(x_bar)
                        y_bar(k) = double(obj.eos(ind(k)).pressure(x_bar(k)));
                    end
                end
            else
                y_bar = NaN;
            end
            
            % Prepare the data: scaled input
            if pr.showscaledinput && ~isempty(obj.eos) && (range(x_tof) > 0)
                x_bar = linspace(min(x_tof), max(x_tof));
                bnorm = obj.betam; % the mass renorm factor
                anorm = obj.alfar; % the radius renorm factor
                if isempty(bnorm), bnorm = nan; end
                if isempty(anorm), anorm = nan; end
                if isscalar(obj.eos)
                    y_bar_scl = double(bnorm/anorm*obj.eos.pressure(x_bar/bnorm));
                else
                    v = 1:length(unique(x_tof));
                    ind = interp1(unique(x_tof), v, x_bar, 'nearest', 'extrap');
                    y_bar_scl = nan(size(x_bar));
                    for k=1:length(x_bar)
                        y_bar_scl(k) = double(...
                            bnorm/anorm*obj.eos(ind(k)).pressure(x_bar(k)/bnorm));
                    end
                end
            else
                y_bar_scl = NaN;
            end
            
            % Plot the lines (pressure in GPa)
            lh(1) = line(x_tof, y_tof/1e9);
            lh(1).LineWidth = 2;
            if isempty(pr.axes)
                lh(1).Color = [0.31, 0.31, 0.31];
            end
            if isempty(obj.name)
                lh(1).DisplayName = 'TOF4 model';
            else
                lh(1).DisplayName = obj.name;
            end
            
            if pr.showinput && any(isfinite(y_bar))
                lh(end+1) = line(x_bar, y_bar/1e9);
                lh(end).Color = 'r';
                lh(end).LineStyle = '--';
                lh(end).DisplayName = 'input barotrope';
            end
            
            if pr.showscaledinput && any(isfinite(y_bar_scl))
                lh(end+1) = line(x_bar, y_bar_scl/1e9);
                lh(end).Color = [0, 0.5, 0];
                lh(end).LineStyle = '--';
                lh(end).DisplayName = 'input barotrope ($\frac{\beta}{\alpha}$-scaled)';
            end
            
            % Style and annotate axes
            if isempty(pr.axes)
                ah.Box = 'on';
                if (range(x_tof) > 0)
                    xlim([min(x_tof),max(x_tof)])
                end
                xlabel('$\rho$ [kg/m$^3$]')
                ylabel('$P$ [GPa]')
            else
                xlim('auto')
            end
            
            % Legend
            legend(ah, 'off')
            gh = legend(ah, 'show','location','nw');
            gh.FontSize = 11;
            
        end
        
        function [ah, lh] = plot_equipotential_surfaces(obj)
            % Visualize a TOFPlanet object by plotting equipotential contours.
            
            % Require R2016a to use the amazing polarplot features
            if verLessThan('matlab','9')
                warning('Equipotential plots require R2016a or later')
                return
            end
            
            % Work on converged objects only
            if isempty(obj.Js)
                warning('There are no equipotentials. Did you run tof.relax_to_HE() yet?')
                return
            end
            
            % Prepare polar axes
            figure;
            ah = polaraxes;
            ah.ThetaZeroLocation = 'top';
            ah.ThetaDir = 'clockwise';
            ah.ThetaAxisUnits = 'rad';
            hold(ah, 'on')
            
            % Plot level surfaces colored by layer density
            cmap = parula;
            rho = obj.rhoi;
            romin = min(rho); romax = max(rho);
            lh = gobjects(size(obj.si));
            for k=1:obj.N
                th = linspace(0,2*pi,60);
                xi = obj.level_surface(obj.si(k)/obj.s0, th);
                lh(k) = polarplot(ah, th, xi);
                lh(k).Tag = 'equisurface';
                if (rho(k) <= romin)
                    ci = 1;
                elseif (rho(k) >= romax)
                    ci = length(cmap);
                else
                    ci = fix((rho(k) - romin)/(romax - romin)*length(cmap)) + 1;
                end
                lh(k).Color = cmap(ci,:);
            end
            
            % Make outer surface more distinct
            lh(1).LineWidth = 2;
            lh(1).Color = 'k';
            
            % Show grid lines above contours
            ah.Layer = 'top';
            
            % Add a colorbar
            ch = colorbar;
            ch.Label.String =...
                sprintf('\\times %.0f kg/m^3', double(max(obj.rhoi)));
            ch.Label.FontSize = 10;
        end
        
    end % End of public methods block
    
    %% Private (or obsolete) methods
    methods (Access = private)
        function beta = renormalize_density(obj)
            beta = obj.mass/obj.M;
            if abs(beta - 1) < 2*eps, return, end % don't renormalize a normal
            obj.betam = beta;
            obj.rhoi = obj.rhoi*beta;
        end
        function y = P_ref(obj)
            U = obj.G*obj.mass/obj.s0^3*obj.si.^2.*obj.Upu();
            rho = 0.5*(obj.rhoi(1:end-1) + obj.rhoi(2:end));
            y = zeros(obj.N, 1)*rho(1)*U(1);
            y(1) = obj.P0;
            y(2:end) = y(1) + cumsum(-rho.*diff(U));
        end
        function y = P_mid(obj)
            % Pressure interpolated to halfway between level surfaces.
            
            v = double(obj.Pi);
            if isempty(v), y = []; return, end
            x = double(obj.si);
            xq = [(x(1:end-1) + x(2:end))/2; x(end)/2];
            y = interp1(x, v, xq, 'pchip')*obj.u.Pa;
        end
        function y = Upu(obj)
            % Following Nettelmann 2017 eqs. B3 and B.4, assuming equipotential.
            s2 = obj.ss.s2;
            s4 = obj.ss.s4;
            A0(obj.N,1) = 0;
            A0 = A0 + obj.SS.S0.*(1 + 2/5*s2.^2 - 4/105*s2.^3 + 2/9*s4.^2 + ...
                43/175*s2.^4 - 4/35*s2.^2.*s4);
            A0 = A0 + obj.SS.S2.*(-3/5*s2 + 12/35*s2.^2 - 234/175*s2.^3 + 24/35*s2.*s4);
            A0 = A0 + obj.SS.S4.*(-5/9*s4 + 6/7*s2.^2);
            A0 = A0 + obj.SS.S0p.*(1);
            A0 = A0 + obj.SS.S2p.*(2/5*s2 + 2/35*s2.^2 + 4/35*s2.*s4 - 2/25*s2.^3);
            A0 = A0 + obj.SS.S4p.*(4/9*s4 + 12/35*s2.^2);
            A0 = A0 + obj.mrot/3*(1 - 2/5*s2 - 9/35*s2.^2 - 4/35*s2.*s4 + 22/525*s2.^3);
            
            y = -A0;
        end
    end % End of private methods block
    
    %% Access methods
    methods
        function val = get.Ui(obj)
            % Following Nettelmann 2017 eqs. B3 and B.4, assuming equipotential.
            if isempty(obj.ss), val = []; return, end
            
            s2 = obj.ss.s2;
            s4 = obj.ss.s4;
            A0(obj.N,1) = 0;
            A0 = A0 + obj.SS.S0.*(1 + 2/5*s2.^2 - 4/105*s2.^3 + 2/9*s4.^2 + ...
                                  43/175*s2.^4 - 4/35*s2.^2.*s4);
            A0 = A0 + obj.SS.S2.*(-3/5*s2 + 12/35*s2.^2 - 234/175*s2.^3 + 24/35*s2.*s4);
            A0 = A0 + obj.SS.S4.*(-5/9*s4 + 6/7*s2.^2);
            A0 = A0 + obj.SS.S0p.*(1);
            A0 = A0 + obj.SS.S2p.*(2/5*s2 + 2/35*s2.^2 + 4/35*s2.*s4 - 2/25*s2.^3);
            A0 = A0 + obj.SS.S4p.*(4/9*s4 + 12/35*s2.^2);
            A0 = A0 + obj.mrot/3*(1 - 2/5*s2 - 9/35*s2.^2 - 4/35*s2.*s4 + 22/525*s2.^3);
            
            val = -obj.G*obj.mass/obj.s0^3*obj.si.^2.*A0;
        end
                
        function val = get.Pi(obj)
            if isempty(obj.Ui) || isempty(obj.rhoi) || isempty(obj.P0)
                val = [];
                return
            end
            U = obj.Ui;
            rho = 0.5*(obj.rhoi(1:end-1) + obj.rhoi(2:end));
            val = zeros(obj.N, 1)*rho(1)*U(1);
            val(1) = obj.P0;
            val(2:end) = val(1) + cumsum(-rho.*diff(U));
        end
        
        function set.eos(obj,val)
            if isempty(val)
                obj.eos = [];
                return
            end
            if ~isa(val,'barotropes.Barotrope')
                error('eos must be a valid instance of class Barotrope')
            end
            assert(numel(val)==1 || numel(val)==obj.N,...
                'eos must be scalar or same length as radius grid.')
            obj.eos = val(:);
        end
        
        function val = get.M(obj)
            if isempty(obj.si) || isempty(obj.rhoi)
                val = [];
            else
                val = obj.total_mass(obj.opts.masmeth);
            end
        end
        
        function val = get.s0(obj)
            if isempty(obj.si)
                val = [];
            else
                val = obj.si(1);
            end
        end
        
        function val = get.rhobar(obj)
            if isempty(obj.M) || isempty(obj.si)
                val = [];
            else
                val = obj.M/(4*pi/3*obj.s0^3);
            end
        end
        
        function val = get.qrot(obj)
            val = obj.mrot*obj.aos^3;
        end
        
        function val = get.a0(obj)
            if isempty(obj.si)
                val = [];
            else
                val = obj.aos*obj.s0;
            end
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
