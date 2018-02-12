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
        bgeos  % optional background barotrope
        fgeos  % optional foreground barotrope
        rhoc   % threshold density indicating transition to solid core
        opts   % holds user configurable options (tip: help tofset)
    end
    properties (SetAccess = private)
        N      % convenience name for length(obj.si)
        ss     % shape functions (returned by tof4.m)
        SS     % shape functions (returned by tof4.m)
        Js     % external gravity coefficients (returned by tof4.m)
        betam  % mass renormalization factor (obj.mass/obj.M)
        alfar  % radius renormalization factor (obj.radius/obj.a0)
    end
    properties (Dependent)
        M      % calculated mass
        mi     % cumulative mass below si
        zi     % mass fraction in heavy elements
        mzi    % cumulative heavy elements mass below si
        M_core % estimated core mass
        R_core % estimated core radius
        M_Z    % estimated total heavy elements mass
        s0     % calculated mean radius (another name for obj.si(1))
        a0     % calculated equatorial radius
        rhobar % calculated mean density
        qrot   % rotation parameter referenced to a0
        NMoI   % normalized moment of inertia
        Ui     % gravitational potential on grid
        Pi     % hydrostatic pressure on grid
        J2     % convenience alias to obj.Js(2)
        J4     % convenience alias to obj.Js(3)
        J6     % convenience alias to obj.Js(4)
        J8     % convenience alias to obj.Js(5)
    end
    properties (GetAccess = private)
        aos    % calculated equatorial to mean radius ratio (from tof4.m)
        G      % Gravitational constant
        u      % let's hold a units struct for convenience
    end
    
    %% A simple constructor
    methods
        function obj = TOFPlanet(varargin)
            % A simple constructor of TOFPlanet objects.
            % TOFPlanet(N, 'OPTION1', VALUE, 'OPTION2', VALUE2,...)
            
            % Populate options struct
            obj.opts = tofset(varargin{:});
            
            % Init privates
            obj.aos = 1;
            try
                if obj.opts.debug
                    obj.u = setUnits;
                else
                    obj.u = setFUnits;
                end
                obj.G = obj.u.gravity;
            catch ME
                if isequal(ME.identifier,'MATLAB:UndefinedFunction')
                    error('Check your path, did you forget to run setws()?\n%s',ME.message)
                end
                rethrow(ME)
            end
        end
    end % End of constructor block
    
    %% Public methods
    methods (Access = public)
        function set_ss_guesses(obj, ss_guesses)
            % Supply a shape functions struct to seed relax_to_HE call.
            
            if nargin < 2, ss_guesses = []; end
            if isempty(ss_guesses), obj.ss = []; return; end
            validateattributes(ss_guesses,{'struct'},{'scalar'},'','ss_guesses',1)
            try
                x = ss_guesses;
                validateattributes(x.s0,{'numeric'},{'column','numel',obj.N},'','s0')
                validateattributes(x.s2,{'numeric'},{'column','numel',obj.N},'','s2')
                validateattributes(x.s4,{'numeric'},{'column','numel',obj.N},'','s4')
                validateattributes(x.s6,{'numeric'},{'column','numel',obj.N},'','s6')
                validateattributes(x.s8,{'numeric'},{'column','numel',obj.N},'','s8')
                obj.ss = ss_guesses;
            catch ME
                warning('ss_guesses not set because:\n%s',ME.message)
            end
        end
        
        function ET = relax_to_barotrope(obj)
            % Relax equilibrium shape functions, Js, and density simultaneously.
            
            % First some checks.
            if isempty(obj.eos)
                warning('TOFPLANET:assertion',...
                    'Set valid barotrope first (<obj>.eos = <barotrope>).')
                return
            end
            if numel(obj.renormalize()) < 2
                warning('TOFPLANET:assertion',...
                    'First set reference mass and equatorial radius.')
                return
            end
            if isempty(obj.mrot)
                warning('TOFPLANET:assertion',...
                    'First set rotation parameter (<obj>.mrot).')
                return
            end
            if isempty(obj.P0)
                warning('TOFPLANET:P0',...
                    'Setting reference pressure to zero (<obj>.P0=0).')
                obj.P0 = 0*obj.u.bar;
            end
            
            % Optional communication
            verb = obj.opts.verbosity;
            if (verb > 0)
                fprintf('Relaxing to desired barotrope...\n\n')
            end
            
            % Ready, set,...
            mihe = obj.opts.MaxIterHE;
            obj.opts.MaxIterHE = 2;
            warning('off','TOF4:maxiter')
            t_rlx = tic;
            
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
            
            % Some clean up
            obj.opts.MaxIterHE = mihe;
            warning('on','TOF4:maxiter')
            
            % Renormalize densities, radii.
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
            
            t_rlx = tic;
            zvec = double(obj.si/obj.si(1));
            dvec = double(obj.rhoi/obj.rhobar);
            if isempty(obj.ss)
                ss_guess = [];
            else
                ss_guess = structfun(@flipud, obj.ss, 'UniformOutput', false);
            end
            [obj.Js, out] = tof4(zvec, dvec, obj.mrot,...
                obj.opts.dJtol, obj.opts.MaxIterHE, ss_guess);
            ET = toc(t_rlx);
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
        
        function ab = renormalize(obj)
            % Match input and calculated mass and equatorial radius.
            try
                a = obj.radius/obj.a0;
                obj.si = obj.si*a;
            catch
                a = [];
            end
            try
                b = obj.mass/obj.M;
                obj.rhoi = obj.rhoi*b;
            catch
                b = [];
            end
            ab = [a, b];
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
        
        function mcore = core_mass(obj, how)
            % Return estimated core mass.
            % mcore = CORE_MASS(how) returns the estimated mass in the planet's
            % core. There are several options for how to decide what if any part
            % of the planet is a core, and these are specified using the argument,
            % how, in the following ways:
            %
            % CORE_MASS('byeos') returns the mass contained in the innermost
            % contiguous block of levels associated with the same Barotrope
            % object.
            %
            % CORE_MASS('bypeaks') uses the peaks function to identify the
            % innermost block of contiguous levels following a detectable density
            % jump.
            %
            % CORE_MASS(rhoc) where x is a numeric scalar returns the mass
            % contained in layers with density greater than rhoc.
            %
            % Note: the TOFPlanet property M_core is obtained with the call
            %       mcore = obj.core_mass(obj.rhoc);
            
            if nargin == 1 && nargout == 0
                help TOFPlanet.core_mass
                return
            end
            if isnumeric(how)
                validateattributes(how,{'numeric','preal'},{'positive','scalar'})
                if how < 0.5*obj.rhobar
                    warning('Are you using the right density units?')
                end
            else
                how = validatestring(how,{'byeos','bypeaks'});
            end
            switch how
                case 'byeos'
                    % Define core as innermost contiguous block of same eos
                    if isempty(obj.eos)
                        error('Object not assigned EOS yet.')
                    end
                    alleos = obj.eos;
                    if isequal(alleos(1), barotropes.ConstDensity(0))
                        zlay = true;
                        alleos(1) = [];
                    else
                        zlay = false;
                    end
                    ind = arrayfun(@isequal, alleos,...
                        repmat(alleos(end), numel(alleos), 1));
                    cind = find(~ind, 1, 'last') + 1;
                    if isempty(cind)
                        mcore = 0;
                    else
                        if zlay, cind = cind + 1; end
                        mcore = obj.mi(cind);
                    end
                case 'bypeaks'
                    % Define core using a density jump
                    dvec = double(obj.rhoi/obj.rhobar)';
                    deltas = [dvec(1), diff(dvec)];
                    cind = peakfinder(deltas);
                    if isempty(obj.rhoc) || obj.rhoi(cind(end)+1) >= obj.rhoc
                        mcore = obj.mi(cind(end));
                    else
                        mcore = 0*obj.M;
                    end
                otherwise
                    % Define core as anything with density over rhoc
                    cind = find(obj.rhoi > how, 1, 'first');
                    mcore = obj.mi(cind);
                    if isempty(mcore), mcore = 0*obj.u.kg; end
            end
        end
        
        function rcore = core_radius(obj, how)
            % Return estimated core radius.
            % rcore = CORE_RADIUS(how) returns the estimated radius of the planet's
            % core. There are several options for how to decide what if any part
            % of the planet is a core, and these are specified using the argument,
            % how, in the following ways:
            %
            % CORE_RADIUS('byeos') returns the radius of the innermost
            % contiguous block of levels associated with the same Barotrope
            % object.
            %
            % CORE_RADIUS('bypeaks') uses the peaks function to identify the
            % innermost block of contiguous levels following a detectable density
            % jump.
            %
            % CORE_radius(rhoc) where x is a numeric scalar returns the radius
            % of the outermost level surface with density greater than rhoc.
            %
            % Note: the TOFPlanet property R_core is obtained with the call
            %       rcore = obj.core_radius(obj.rhoc);
            
            if nargin == 1 && nargout == 0
                help TOFPlanet.core_radius
                return
            end
            if isnumeric(how)
                validateattributes(how,{'numeric','preal'},{'positive','scalar'})
                if how < 0.5*obj.rhobar
                    warning('Are you using the right density units?')
                end
            else
                how = validatestring(how,{'byeos','bypeaks'});
            end
            switch how
                case 'byeos'
                    % Define core as innermost contiguous block of same eos
                    if isempty(obj.eos)
                        error('Object not assigned EOS yet.')
                    end
                    alleos = obj.eos;
                    if isequal(alleos(1), barotropes.ConstDensity(0))
                        zlay = true;
                        alleos(1) = [];
                    else
                        zlay = false;
                    end
                    ind = arrayfun(@isequal, alleos,...
                        repmat(alleos(end), numel(alleos), 1));
                    cind = find(~ind, 1, 'last') + 1;
                    if isempty(cind)
                        rcore = 0;
                    else
                        if zlay, cind = cind + 1; end
                        rcore = obj.si(cind);
                    end
                case 'bypeaks'
                    % Define core using a density jump
                    dvec = double(obj.rhoi/obj.rhobar)';
                    deltas = [dvec(1), diff(dvec)];
                    cind = peakfinder(deltas);
                    rcore = obj.si(cind(end));
                otherwise
                    cind = find(obj.rhoi > how, 1, 'first');
                    rcore = obj.si(cind);
                    if isempty(rcore), rcore = 0*obj.u.m; end
            end
        end
        
        function [ah, lh, gh] = plot_barotrope(obj, varargin)
            % Plot P(rho) of current model and optionally of input barotrope.
            
            % Don't bother if there is no pressure
            if isempty(obj.Pi)
                warning('Uninitialized object. Remember to set obj.P0?')
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
        
        function [ah, lh, gh] = plot_rho_of_r(obj, varargin)
            % Plot rho(r) where r is mean radius.
            
            % Don't bother if uninitialized
            if isempty(obj.rhobar)
                warning('Uninitialized object.')
                return
            end
            
            % Input parsing
            p = inputParser;
            p.FunctionName = mfilename;
            p.addParameter('axes', [], @(x)isscalar(x) && isgraphics(x, 'axes'))
            p.addParameter('plottype', 'line', @(x)isrow(x) && ischar(x))
            p.addParameter('removeeos',[],@(x)isscalar(x)&&isa(x,'barotropes.Barotrope'))
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
            
            % Prepare the data
            x = [double(obj.si/obj.s0); 0];
            y = double([obj.rhoi; obj.rhoi(end)]);
            if ~isempty(pr.removeeos) % optionally remove background density
                if isempty(obj.Pi)
                    warning('Uninitialized object. Remember to set obj.P0?')
                    return
                end
                bgrho = pr.removeeos.density(double(obj.Pi));
                bgrho(isnan(bgrho)) = 0;
                bgrho = [bgrho; bgrho(end)];
                y = max(y - bgrho, 0);
            end
            
            % Plot the lines (density in 1000 kg/m^3)
            if isequal(lower(pr.plottype), 'stairs')
                lh = stairs(x, y/1000);
            elseif isequal(lower(pr.plottype), 'line')
                lh = line(x, y/1000);
            else
                error('Unknown value of parameter plottype.')
            end
            lh.LineWidth = 2;
            if isempty(pr.axes)
                lh.Color = [0.31, 0.31, 0.31];
            end
            if isempty(obj.name)
                lh.DisplayName = 'TOF model';
            else
                lh.DisplayName = obj.name;
            end
            
            % Style and annotate axes
            if isempty(pr.axes)
                ah.Box = 'on';
                xlabel('Level surface radius, $s/s_0$', 'fontsize', 12)
                ylabel('$\rho$ [1000 kg/m$^3$]', 'fontsize', 12)
                if ~isempty(pr.removeeos)
                    s = sprintf('$\\rho - \\rho_{\\mathrm{%s}}$ [1000 kg/m$^3$]',pr.removeeos.name);
                    ylabel(s, 'fontsize', 12)
                end
            else
                xlim('auto')
            end
            
            % Legend
            legend(ah, 'off')
            gh = legend(ah, 'show','location','ne');
            gh.FontSize = 11;
        end
        
        function [ah, lh, gh] = plot_P_of_r(obj, varargin)
            % Plot P(r) where r is mean radius.
            
            % Don't bother if uninitialized
            if isempty(obj.rhobar)
                warning('Uninitialized object.')
                return
            end
            
            % Input parsing
            p = inputParser;
            p.FunctionName = mfilename;
            p.addParameter('axes', [], @(x)isscalar(x) && isgraphics(x, 'axes'));
            p.addParameter('pressurepoint', 'top', @(x)isrow(x) && ischar(x));
            p.addParameter('plottype', 'stairs', @(x)isrow(x) && ischar(x));
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
            
            % Prepare the data
            x = double(obj.si/obj.s0);
            P = double(obj.Pi);
            P_c = interp1(x, P, 0, 'pchip');
            x = [x; 0];
            y = [P; P_c];
            
            % Plot the lines (pressure in GPa)
            if isequal(lower(pr.plottype), 'stairs')
                lh = stairs(x, y/1e9);
            elseif isequal(lower(pr.plottype), 'line')
                lh = line(x, y/1e9);
            else
                error('Unknown value of parameter plottype.')
            end
            lh.LineWidth = 2;
            if isempty(pr.axes)
                lh.Color = [0.31, 0.31, 0.31];
            end
            if isempty(obj.name)
                lh.DisplayName = 'TOF model';
            else
                lh.DisplayName = obj.name;
            end
            
            % Style and annotate axes
            if isempty(pr.axes)
                ah.Box = 'on';
                xlabel('Level surface radius, $s/s_0$', 'fontsize', 12)
                ylabel('$P$ [GPa]', 'fontsize', 12)
            else
                xlim('auto')
            end
            
            % Legend
            legend(ah, 'off')
            gh = legend(ah, 'show','location','ne');
            gh.FontSize = 11;
        end
        
        function T = report_card(obj, obs)
            % REPORT_CARD Table summary of model's vital statistics.
            
            % Minimal checks
            narginchk(1,2);
            try
                obj.J2;
            catch
                warning('Uncooked object.')
                return
            end
            
            % Basic table
            vitals = {'Mass [kg]'; 'J2'; 'J4'; 'J6'; 'J8'; 'NMoI'; '"core" mass [kg]'};
            TOF1 = [obj.M; obj.J2; obj.J4; obj.J6; obj.J8; obj.NMoI; obj.M_core];
            TOF1 = double(TOF1);
            if isempty(obj.M_core), TOF1 = [TOF1; NaN]; end
            T = table(TOF1, 'RowNames', vitals);
            if ~isempty(obj.name)
                vname = matlab.lang.makeValidName(obj.name);
                T.Properties.VariableNames{'TOF1'} = vname;
            end
            if nargin == 1, return, end
            
            % Optionally compare with something
            try
                oM = obs.M;
            catch
                oM = NaN;
            end
            try
                oJ2 = obs.J2;
                oJ4 = obs.J4;
                oJ6 = obs.J6;
                oJ8 = obs.J8;
            catch
                oJ2 = NaN;
                oJ4 = NaN;
                oJ6 = NaN;
                oJ8 = NaN;
            end
            try
                oNMoI = obs.NMoI;
            catch
                oNMoI = NaN;
            end
            try
                oM_core = obs.M_core;
            catch
                oM_core = NaN;
            end
            try
                oname = obs.name;
            catch
                oname = [];
            end
            OBS1 = [oM; oJ2; oJ4; oJ6; oJ8; oNMoI; oM_core];
            OBS1 = double(OBS1);
            if isempty(oM_core), OBS1 = [OBS1; NaN]; end
            T = [T table(OBS1)];
            if ~isempty(oname)
                vname = matlab.lang.makeValidName(obs.name);
                try
                    T.Properties.VariableNames{'OBS1'} = vname;
                catch
                    T.Properties.VariableNames{'OBS1'} = ['x_',vname];
                end
            end
            DIFF = (TOF1 - OBS1)./TOF1;
            T = [T table(DIFF, 'VariableNames', {'diff'})];
        end
        
        function s = to_struct(obj, rdc, keepss)
            % Convert object to static struct keeping only essential fields.
            
            if nargin < 2, rdc = 1; end % 0=none, 1=to double, 2=to single, 3=to scalars
            if nargin < 3, keepss = false; end % keep ss e.g. to help essample
            
            s.name   = obj.name;
            s.M      = obj.M;
            s.s0     = obj.s0;
            s.a0     = obj.a0;
            s.rhoc   = obj.rhoc;
            s.M_core = obj.M_core;
            s.R_core = obj.R_core;
            s.M_Z    = obj.M_Z;
            s.rhobar = obj.rhobar;
            s.mrot   = obj.mrot;
            s.qrot   = obj.qrot;
            s.J2     = obj.J2;
            s.J4     = obj.J4;
            s.J6     = obj.J6;
            s.J8     = obj.J8;
            s.NMoI   = obj.NMoI;
            s.si     = obj.si;
            s.rhoi   = obj.rhoi;
            s.Pi     = obj.Pi;
            s.mi     = obj.mi;
            s.zi     = obj.zi;
            
            if rdc == 1
                s = structfun(@double, s, 'UniformOutput', false);
                s.name = obj.name;
            end
            if rdc == 2
                s = structfun(@single, s, 'UniformOutput', false);
                s.name = obj.name;
            end
            if rdc == 3
                s = structfun(@single, s, 'UniformOutput', false);
                s.name = obj.name;
                s.si     = [];
                s.rhoi   = [];
                s.Pi     = [];
                s.mi     = [];
                s.zi     = [];
            end
            
            try
                s.eos = obj.eos.name;
            catch
                s.eos = '';
            end
            try
                s.bgeos = obj.bgeos.name;
            catch
                s.bgeos = '';
            end
            try
                s.fgeos = obj.fgeos.name;
            catch
                s.fgeos = '';
            end
            
            if keepss
                s.ss = obj.ss;
            end
        end
        
        function T = to_table(obj)
            % Create table of the grid quantities.
            
            T = table;
            T.si = double(obj.si);
            T.rhoi = double(obj.rhoi);
            T.Pi = double(obj.Pi);
            T.mi = double(obj.mi);
            if ~isempty(obj.zi)
                T.zi = double(obj.zi);
            end
        end
        
        function to_ascii(obj, fname)
            % Export the state of the model as ascii file.
            
            % File name
            if nargin == 1, fname = obj.name; end
            if isempty(fname), fname = 'model1.txt'; end
            validateattributes(fname, {'char'}, {'row'}, '', 'fname', 1)
            
            % Open file
            fid = fopen(fname,'wt');
            cleanup = onCleanup(@()fclose(fid));
            
            % Write the header
            fprintf(fid,'# Rotating fluid planet modeled 4th-order Theory of Figures.\n');
            fprintf(fid,'#\n');
            fprintf(fid,'# Model name: %s\n', obj.name);
            fprintf(fid,'#\n');
            fprintf(fid,'# Scalar quantities:\n');
            fprintf(fid,'# N layers = %d\n',obj.N);
            fprintf(fid,'# Mass M = %g kg\n', double(obj.M));
            fprintf(fid,'# Mean radius       s0 = %0.6e m\n', double(obj.s0));
            fprintf(fid,'# Equatorial radius a0 = %0.6e m\n', double(obj.a0));
            fprintf(fid,'# Rotation parameter m = %0.6f\n', double(obj.mrot));
            fprintf(fid,'# Rotation parameter q = %0.6f\n', double(obj.qrot));
            fprintf(fid,'# Core mass fraction M_core/M = %g\n', ...
                double(obj.M_core)/double(obj.M));
            fprintf(fid,'#\n');
            fprintf(fid,'# Calculated gravity zonal harmonics (x 10^6):\n');
            fprintf(fid,'# J2  = %12.6f\n', obj.J2*1e6);
            fprintf(fid,'# J4  = %12.6f\n', obj.J4*1e6);
            fprintf(fid,'# J6  = %12.6f\n', obj.J6*1e6);
            fprintf(fid,'# J8  = %12.6f\n', obj.J8*1e6);
            fprintf(fid,'#\n');
            fprintf(fid,'# Column data description (MKS):\n');
            fprintf(fid,'# i     - level surface index (increasing with depth)\n');
            fprintf(fid,'# s_i   - mean radius of level surface i\n');
            fprintf(fid,'# rho_i - density on level surfaces i\n');
            fprintf(fid,'# P_i   - pressure on level surface i\n');
            fprintf(fid,'# m_i   - mass below level surface i\n');
            fprintf(fid,'#\n');
            
            % Write the data
            fprintf(fid,'# Column data:\n');
            fprintf(fid,'# %-4s  ','i');
            fprintf(fid,'%-10s  ','s_i');
            fprintf(fid,'%-7s  ','rho_i');
            fprintf(fid,'%-10s  ','P_i','m_i');
            fprintf(fid,'\n');
            for k=1:obj.N
                fprintf(fid,'  %-4d  ',k);
                fprintf(fid,'%10.4e  ', double(obj.si(k)));
                fprintf(fid,'%7.1f  ', double(obj.rhoi(k)));
                fprintf(fid,'%10.4e  ', double(obj.Pi(k)), double(obj.mi(k)));
                fprintf(fid,'\n');
            end
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
        function set.name(obj,val)
            if ~isempty(val)
                validateattributes(val, {'char'}, {'row'})
            end
            obj.name = val;
        end
        
        function set.si(obj, val)
            assert(isnumeric(val) && isvector(val) && ~any(val<0),...
                'obj.si must be a nonnegative vector.')
            assert(all(diff(val)<=0),'obj.si must be non-ascending.')
            obj.si = val(:);
        end
        
        function set.rhoi(obj, val)
            assert(isnumeric(val) && isvector(val),...
                'obj.rhoi must be a nonnegative vector.')
            if any(val<0)
                warning('TOFPLANET:assertion','negative density. Is this on purpose?')
            end
            obj.rhoi = val(:);
        end
        
        function val = get.Ui(obj)
            % Following Nettelmann 2017 eqs. B3 and B.4, assuming equipotential.
            if isempty(obj.ss), val = []; return, end
            if isempty(obj.SS), val = []; return, end
            
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
            obj.eos = val(:);
        end
        
        function set.bgeos(obj,val)
            if isempty(val)
                obj.bgeos = [];
                return
            end
            if ~isa(val,'barotropes.Barotrope')
                error('bgeos must be a valid instance of class Barotrope')
            end
            obj.bgeos = val;
        end

        function set.fgeos(obj,val)
            if isempty(val)
                obj.fgeos = [];
                return
            end
            if ~isa(val,'barotropes.Barotrope')
                error('fgeos must be a valid instance of class Barotrope')
            end
            obj.fgeos = val;
        end
        
        function val = get.M(obj)
            if isempty(obj.si) || isempty(obj.rhoi)
                val = [];
            else
                val = obj.total_mass(obj.opts.masmeth);
            end
        end
        
        function set.mass(obj,val)
            validateattributes(val,{'numeric'},{'positive','scalar'})
            obj.mass = val;
        end
        
        function set.radius(obj,val)
            validateattributes(val,{'numeric'},{'positive','scalar'})
            obj.radius = val;
        end
        
        function val = get.mi(obj)
            % mass _below_ level i
            if isempty(obj.si) || isempty(obj.rhoi)
                val = [];
            else
                rho = obj.rhoi;
                s = obj.si;
                val(obj.N) = 4*pi/3*rho(obj.N)*s(obj.N)^3;
                for k=obj.N-1:-1:1
                    val(k) = val(k+1) + 4*pi/3*rho(k)*(s(k)^3 - s(k+1)^3);
                end
                val = val';
            end
        end
        
        function val = get.mzi(obj)
            % heavy element mass below level i
            if isempty(obj.si) || isempty(obj.rhoi) || isempty(obj.zi)
                val = [];
            else
                rho = obj.rhoi;
                s = obj.si;
                z = obj.zi;
                val(obj.N) = 4*pi/3*rho(obj.N)*s(obj.N)^3*z(obj.N);
                for k=obj.N-1:-1:1
                    cz = max(min(z(k), 1), 0);
                    val(k) = val(k+1) + 4*pi/3*rho(k)*(s(k)^3 - s(k+1)^3)*cz;
                end
                val = val';
            end
        end
        
        function val = get.zi(obj)
            % heavy element mass fraction on level i
            if isempty(obj.bgeos) || isempty(obj.fgeos) || isempty(obj.Pi)
                val = [];
            else
                val = nan(obj.N,1);
                for k=1:obj.N
                    roxy = obj.bgeos.density(obj.Pi(k));
                    roz = obj.fgeos.density(obj.Pi(k));
                    ro = obj.rhoi(k);
                    val(k) = (1/ro - 1/roxy)/(1/roz - 1/roxy);
                    if isnan(val(k)), val(k) = 0; end
                end
            end
        end
        
        function val = get.M_Z(obj)
            try
                val = obj.mzi(1);
            catch
                val = [];
            end
        end
        
        function val = get.M_core(obj)
            try
                val = obj.core_mass(obj.rhoc);
            catch
                val = [];
            end
        end
        
        function val = get.R_core(obj)
            try
                val = obj.core_radius(obj.rhoc);
            catch
                val = [];
            end
        end
        
        function val = get.s0(obj)
            if isempty(obj.si)
                val = [];
            else
                val = obj.si(1);
            end
        end
        
        function val = get.N(obj)
            if isempty(obj.si) || isempty(obj.rhoi)
                val = 0;
            elseif length(obj.si) == length(obj.rhoi)
                val = length(obj.si);
            else
                error('length(si) = %g ~= length(rhoi) = %g',...
                    length(obj.si),length(obj.rhoi))
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
        
        function val = get.NMoI(obj)
            if isempty(obj.rhobar)
                val = [];
            else
                rho = obj.rhoi;
                s = obj.si;
                for k=1:obj.N-1
                    m = 4*pi/3*rho(k)*(s(k)^3 - s(k+1)^3);
                    I(k) = 2/5*m*(s(k)^5 - s(k+1)^5)/(s(k)^3 - s(k+1)^3); %#ok<AGROW>
                end
                val = sum(I)/(obj.M*obj.s0^2);
            end
        end
        
        function val = get.J2(obj)
            if isempty(obj.Js)
                val = 0;
            else
                val = obj.Js(2);
            end
        end
        
        function val = get.J4(obj)
            if isempty(obj.Js)
                val = 0;
            else
                val = obj.Js(3);
            end
        end
        
        function val = get.J6(obj)
            if isempty(obj.Js)
                val = 0;
            else
                val = obj.Js(4);
            end
        end
        
        function val = get.J8(obj)
            if isempty(obj.Js)
                val = 0;
            else
                val = obj.Js(5);
            end
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
