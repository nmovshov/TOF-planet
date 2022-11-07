classdef TOFPlanet < handle
    %TOFPLANET Interior model of rotating fluid planet.
    %   This class implements a model of a rotating fluid planet using Theory of
    %   Figures to calculate the hydrostatic equilibrium shape and resulting
    %   gravity field. A TOFPlanet object is defined by a densiy profile rho(s),
    %   supplied by the user and stored in the column vectors obj.si and
    %   obj.rhoi, indexed from the surface in. To complete the definition the
    %   user must also specify a mass, equatorial radius, and rotation period.
    %   With these a gravity field and equilibrium shape can be determined,
    %   with a call to obj.relax_to_HE(). Note, however, that the oblate shape
    %   calculated with relax_to_HE() preserves the mass and mean radius of the
    %   planet, but not the equatorial radius. A call to fix_radius()
    %   renormalizs the si vector to match the reference equatorial radius, at
    %   the cost of modifying the implied mass. A call to renormalize()
    %   modifies both si and rhoi to preserve the reference mass and equatorial
    %   radius, at the cost of modifying the assigned density. It is not
    %   possible to define mass, radius, and density simultaneously.
    %
    %   Alternatively the user may supply a barotrope object, stored in
    %   obj.eos, and call obj.relax_to_barotrope() to iteratively find a
    %   density profile consistent with the calculated equilibrium pressure. To
    %   this end a boundary pressure must also be given, in obj.P0, and an
    %   initial density profile guess is still required (can be a simple one).
    %   Again, we can't simultaneously impose an exact mass, radius, and
    %   barotrope. By default the reference mass and radius will be honored, by
    %   renormalizing the converged density profile, and this will modify the
    %   _effective_ barotrope (to override set obj.opts.renorm=false).
    
    %% Properties
    properties
        name   % model name
        mass   % reference mass
        radius % reference radius (equatorial!)
        period % reference rotation period
        P0     % reference pressure
        si     % vector of mean radii (top down, s0=si(1) is outer radius)
        rhoi   % vector of densities on si grid
        eos    % barotrope(s) (tip: help barotropes for options)
        bgeos  % optional background barotrope
        fgeos  % optional foreground barotrope
        opts   % holds user configurable options (tip: help tofset)
    end
    properties (SetAccess = private)
        N      % convenience name for length(obj.si)
        ss     % shape functions (returned by tof<n>.m)
        SS     % shape functions (returned by tof<n>.m)
        A0     % dimensionless potential (returned by tof<n>.m)
        Js     % external gravity coefficients (returned by tof<n>.m)
        betam  % mass renormalization factor returned by obj.renormalize()
        alfar  % radius renormalization factor returned by obj.renormalize()
    end
    properties (Dependent)
        M      % calculated mass (equal to mass *after* renorm)
        a0     % calculated equatorial radius (equal to radius *after* renorm)
        mi     % cumulative mass below si
        ai     % equatorial radii on level surfaces
        s0     % surface mean radius (another name for obj.si(1))
        rhobar % calculated mean density
        wrot   % rotation frequency, 2pi/period
        qrot   % rotation parameter wrot^2a0^3/GM
        mrot   % rotation parameter, wrot^2s0^3/GM
        J2     % convenience alias to obj.Js(2)
        J4     % convenience alias to obj.Js(3)
        J6     % convenience alias to obj.Js(4)
        J8     % convenience alias to obj.Js(5)
        J10    % convenience alias to obj.Js(6)
        J12    % convenience alias to obj.Js(7)
    end
    properties (GetAccess = private)
        aos    % calculated equatorial to mean radius ratio (from tof<n>.m)
        G      % Gravitational constant
    end
    
    %% A simple constructor
    methods
        function obj = TOFPlanet(varargin)
            % A simple constructor of TOFPlanet objects.
            % TOFPlanet(N, OBS, 'OPTION1', VALUE, 'OPTION2', VALUE2,...)
            
            % If first argument is observables struct use it to set references
            if isa(varargin{1}, 'observables')
                obj.set_observables(varargin{1});
                varargin(1) = [];
            end
            % Populate options struct
            obj.opts = tofset(varargin{:});
            
            % Init privates
            obj.aos = 1;
            obj.G = 6.67430e-11; % m^3 kg^-1 s^-2 (2018 NIST reference)
        end
    end % End of constructor block
    
    %% Public methods
    methods (Access = public)
        function obj = set_ss_guesses(obj, ss_guesses)
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
                warning(ME.identifier,'ss_guesses not set because:\n%s',ME.message)
            end
        end
        
        function obj = set_observables(obj, obs)
            % Copy physical properties from an +observables struct.
            obj.mass = obs.M;
            obj.radius = obs.a0;
            obj.period = obs.P;
            obj.P0 = obs.P0;
            try
                obj.bgeos = obs.bgeos;
                obj.fgeos = obs.fgeos;
            catch
            end
        end

        function ET = relax_to_rotation(obj)
            % Relax equilibrium shape and rotation period simultaneously.
            
            % First some checks.
            if isempty(obj.si) || isempty(obj.rhoi)
                warning('TOFPLANET:assertion',...
                    'First radius and density vectors (<obj>.si,<obj>.rhoi).')
                return
            end
            if numel(obj.renormalize()) < 2
                warning('TOFPLANET:assertion',...
                    'First set reference mass and equatorial radius.')
                return
            end
            if isempty(obj.mrot)
                warning('TOFPLANET:assertion',...
                    'First set rotation period (<obj>.period).')
                return
            end
            
            % Optional communication
            verb = obj.opts.verbosity;
            if (verb > 0)
                fprintf('Relaxing to reference rotation period...\n\n')
            end
            
            % Ready, set,...
            warning('off','TOF4:maxiter')
            warning('off','TOF7:maxiter')
            if obj.opts.toforder == 4
                tofun = @tof4;
            else
                tofun = @tof7;
            end
            t_rlx = tic;
            
            % Main loop
            iter = 1;
            while (iter <= obj.opts.MaxIterRot)
                t_pass = tic;
                if (verb > 0)
                    fprintf('Rotationpass %d (of max %d)...\n',...
                        iter, obj.opts.MaxIterRot)
                end
                
                old_Js = obj.Js;
                if isempty(old_Js)
                    old_Js = [-1, zeros(1,obj.opts.toforder)];
                end
                old_m = obj.mrot;
                
                % Call the tof algorithm
                if isempty(obj.ss)
                    ss_guesses = struct();
                else
                    ss_guesses = structfun(@flipud, obj.ss, 'UniformOutput', false);
                end
                [obj.Js, out] = tofun(obj.si, obj.rhoi, obj.mrot,...
                    'tol',obj.opts.dJtol, 'maxiter',obj.opts.MaxIterHE,...
                    'xlevels',obj.opts.xlevels, 'ss_guesses',ss_guesses);
                obj.ss = structfun(@flipud, out.ss, 'UniformOutput', false);
                obj.SS = structfun(@flipud, out.SS, 'UniformOutput', false);
                obj.aos = out.a0;
                obj.A0 = flipud(out.A0);
                
                % Renormalize to reference mass and radius
                obj.renormalize();
                
                % Calculate changes in shape/density
                dJs = abs((obj.Js - old_Js)./old_Js);
                dJs = max(dJs(isfinite(dJs)));
                drot = abs(obj.mrot - old_m);
                
                if (verb > 0)
                    fprintf('Rotationpass %d (of max %d)...done. (%g sec)\n',...
                        iter, obj.opts.MaxIterRot, toc(t_pass))
                    fprintf('var(drot) = %g (%g required); dJ = %g (%g required).\n\n',...
                        drot, obj.opts.drottol, dJs, obj.opts.dJtol)
                end
                
                % The stopping criterion is to satisfy both J and rho tolerance
                if (drot <= obj.opts.drottol) && dJs < obj.opts.dJtol
                    break
                end
                
                % end the main loop
                iter = iter + 1;
            end
            ET = toc(t_rlx);
            if iter > obj.opts.MaxIterRot
                warning('TOFPLANET:maxiter','Rotation may not be fully converged.')
            end
            
            % Renorm and record factors
            % TODO: if we still need betam we must fix this
            renorms = obj.renormalize();
            obj.alfar = renorms(1);
            obj.betam = renorms(2);
            
            % Some clean up
            warning('on','TOF4:maxiter')
            warning('on','TOF7:maxiter')
            
            % Optional communication
            if (verb > 0)
                fprintf('Relaxing to reference rotation...done.\n')
                fprintf('Total elapsed time %s\n',lower(seconds2human(ET)))
            end
        end

        function ET = relax_to_barotrope(obj)
            % Relax equilibrium shape, gravity, and density simultaneously.
            
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
                    'First set rotation period (<obj>.period).')
                return
            end
            if isempty(obj.P0)
                warning('TOFPLANET:P0',...
                    'First set reference pressure (<obj>.P0).')
                return
            end
            
            % Optional communication
            verb = obj.opts.verbosity;
            if (verb > 0)
                fprintf('Relaxing to desired barotrope...\n\n')
            end
            
            % Ready, set,...
            warning('off','TOF4:maxiter')
            warning('off','TOF7:maxiter')
            if obj.opts.toforder == 4
                tofun = @tof4;
            else
                tofun = @tof7;
            end
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
                if isempty(old_Js)
                    old_Js = [-1, zeros(1,obj.opts.toforder)];
                end
                old_ro = obj.rhoi;
                
                % Call the tof algorithm
                if isempty(obj.ss)
                    ss_guesses = struct();
                else
                    ss_guesses = structfun(@flipud, obj.ss, 'UniformOutput', false);
                end
                [obj.Js, out] = tofun(obj.si, obj.rhoi, obj.mrot,...
                    'tol',obj.opts.dJtol, 'maxiter',obj.opts.MaxIterHE,...
                    'xlevels',obj.opts.xlevels, 'ss_guesses',ss_guesses);
                obj.ss = structfun(@flipud, out.ss, 'UniformOutput', false);
                obj.SS = structfun(@flipud, out.SS, 'UniformOutput', false);
                obj.aos = out.a0;
                obj.A0 = flipud(out.A0);
                
                % Update density from barotrope and renormalize
                obj.update_densities();
                obj.renormalize();
                
                % Calculate changes in shape/density
                dJs = abs((obj.Js - old_Js)./old_Js);
                dJs = max(dJs(isfinite(dJs)));
                dro = obj.rhoi./old_ro;
                dro = var(dro(isfinite(dro)));
                
                if (verb > 0)
                    fprintf('Baropass %d (of max %d)...done. (%g sec)\n',...
                        iter, obj.opts.MaxIterBar, toc(t_pass))
                    fprintf('var(drho) = %g (%g required); dJ = %g (%g required).\n\n',...
                        dro, obj.opts.drhotol, dJs, obj.opts.dJtol)
                end
                
                % The stopping criterion is to satisfy both J and rho tolerance
                if (dro < obj.opts.drhotol) && dJs < obj.opts.dJtol
                    break
                end
                
                % end the main loop
                iter = iter + 1;
            end
            ET = toc(t_rlx);
            if iter > obj.opts.MaxIterBar
                warning('TOFPLANET:maxiter','Pressure/density may not be fully converged.')
            end
            
            % Renorm and record factors
            % TODO: if we still need betam we must fix this
            renorms = obj.renormalize();
            obj.alfar = renorms(1);
            obj.betam = renorms(2);
            
            % Some clean up
            warning('on','TOF4:maxiter')
            warning('on','TOF7:maxiter')
            
            % Optional communication
            if (verb > 0)
                fprintf('Relaxing to desired barotrope...done.\n')
                fprintf('Total elapsed time %s\n',lower(seconds2human(ET)))
            end
        end
        
        function [ET, dJ] = relax_to_HE(obj)
            % Call tof<n> once to obtain equilibrium shape and gravity.
            
            if (obj.opts.verbosity > 1)
                fprintf('  Relaxing to hydrostatic equilibrium...\n')
            end
            
            t_rlx = tic;
            zvec = obj.si/obj.si(1);
            dvec = obj.rhoi/obj.rhobar;
            if isempty(obj.ss)
                ss_guess = struct();
            else
                ss_guess = structfun(@flipud, obj.ss, 'UniformOutput', false);
            end
            if obj.opts.toforder == 4
                tofun = @tof4;
            else
                tofun = @tof7;
            end
            [obj.Js, out] = tofun(zvec, dvec, obj.mrot,...
                'tol',obj.opts.dJtol, 'maxiter',obj.opts.MaxIterHE,...
                'ss_guesses',ss_guess, 'xlevels',obj.opts.xlevels);
            ET = toc(t_rlx);
            dJ = out.dJs;
            obj.ss = structfun(@flipud, out.ss, 'UniformOutput', false);
            obj.SS = structfun(@flipud, out.SS, 'UniformOutput', false);
            obj.aos = out.a0;
            obj.A0 = flipud(out.A0);
            
            if (obj.opts.verbosity > 1)
                fprintf('  Relaxing to hydrostatic equilibrium...done.\n')
                fprintf('  Elapsed time %s\n',lower(seconds2human(ET)))
            end
        end
        
        function dro = update_densities(obj)
            % Set level surface densities to match prescribed barotrope.
            
            if isempty(obj.eos)
                error('TOFPLANET:noeos','Make sure input barotrope (<obj>.eos) is set.')
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
        
        function obj = fix_radius(obj)
            % Resize planet to match equatorial radius to observed value.
            
            if isempty(obj.radius) || isempty(obj.a0) || isempty(obj.si)
                warning('TOFPLANET:noref','Missing information; no action.')
                return
            end
            obj.si = obj.si*obj.radius/obj.a0;
        end
        
        function obj = fix_mass(obj)
            % Rescale density to match planet mass to observed value.
            
            if isempty(obj.mass) || isempty(obj.rhoi)
                warning('TOFPLANET:noref','Missing information; no action.')
                return
            end
            obj.rhoi = obj.rhoi*obj.mass/obj.M;
        end
        
        function r = level_surface(obj, el, theta)
            % Return r_l(theta) by expansion of Legendre polynomials.
            
            validateattributes(el,{'numeric'},{'nonnegative','scalar','<=',1},'','l',1)
            validateattributes(theta,{'numeric'},{'vector','>=',0,'<=',2*pi},'','theta',2)
            
            k = find(obj.si/obj.si(1) <= el, 1, 'first');
            mu = cos(theta);
            shp = obj.ss.s0(k)*Pn(0,mu) + obj.ss.s2(k)*Pn(2,mu) + ...
                obj.ss.s4(k)*Pn(4,mu) + obj.ss.s6(k)*Pn(6,mu) + ...
                obj.ss.s8(k)*Pn(8,mu);
            r = obj.si(k)*(1 + shp);
        end
        
        function s = verify_mean_radii(obj)
            % Sanity check: relative error of calculated vs. assigned si.
            %
            % Volume-integrate each level surfaces to check that the mean
            % radius comes out right. If the ToF is implemented correctly they
            % should be very close.
            
            s = nan(obj.N,1);
            for k=1:obj.N
                rmax = @(teta)obj.level_surface(obj.si(k)/obj.si(1), teta);
                fun = @(teta)(rmax(teta).^3).*sin(teta);
                V = 2*pi/3*integral(fun, 0, pi);
                s(k) = ((V/(4*pi/3))^(1/3) - obj.si(k))/obj.si(k);
            end
        end
        
        function I = NMoI(obj, reduce)
            % Return moment of inertia normalized by a0.
            
            if nargin < 2 || isempty(reduce), reduce = 'sum'; end
            reduce = validatestring(reduce, {'sum', 'csum', 'none'});
            
            if strcmp(obj.opts.moimeth, 'layerz')
                deltas = [obj.rhoi(1); diff(obj.rhoi)];
            elseif strcmp(obj.opts.moimeth, 'midlayerz')
                romid = [(obj.rhoi(1:end-1) + obj.rhoi(2:end))/2; obj.rhoi(end)];
                deltas = [romid(1); diff(romid)];
            else
                error('Unknown moimeth %s',obj.opts.moimeth) %#ok<ERTAG> 
            end
            num(obj.N) = 0;
            den = 0;
            [mus, gws] = gauleg(0, 1, 48); % Abscissas and weights for Gauss-quad
            tetas = acos(mus);
            p2term = 1 - Pn(2, mus);
            for k=1:obj.N
                if isempty(obj.ss)
                    ximax = obj.si(k)/obj.a0*ones(size(tetas));
                else
                    ximax = obj.level_surface(obj.si(k)/obj.si(1),tetas)/obj.a0;
                end
                fun1 = deltas(k)*(ximax.^5).*p2term;
                fun2 = deltas(k)*(ximax.^3);
                num(k) = fun1*gws';
                den = den + fun2*gws';
            end
            if isequal(reduce, 'none'), I = (2/5)*(num/den); end
            if isequal(reduce, 'sum'), I = (2/5)*sum(num)/den; end
            if isequal(reduce, 'csum'), I = (2/5)*cumsum(num)/den; end
        end
        
        function m = total_mass(obj,method)
            if nargin == 1, method = 'layerz'; end
            switch lower(method)
                case 'trapz'
                    m = -4*pi*trapz(double(obj.si), double(obj.rhoi.*obj.si.^2));
                case 'layerz'
                    drho = [obj.rhoi(1); diff(obj.rhoi)];
                    m = (4*pi/3)*sum(drho.*obj.si.^3);
                case 'integralz'
                    x = double(obj.si);
                    v = double(obj.rhoi.*obj.si.^2);
                    fun = @(z)interp1(x, v, z, 'pchip');
                    m = 4*pi*integral(fun, 0 , x(1));
                otherwise
                    error('Unknown mass calculation method.') %#ok<ERTAG> 
            end
        end
        
    end % End of public methods block
    
    %% Visualizers
    methods (Access = public)
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
                lh(1).DisplayName = 'TOF model';
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
                error('Unknown value of parameter plottype.') %#ok<ERTAG> 
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
        
        function [ah, lh, gh] = plot_residual_rho_of_r(obj, varargin)
            % Plot rho(r)-rho_xy(r) using a background eos.
            
            % Don't bother if uninitialized
            if isempty(obj.rhobar)
                warning('Uninitialized object.')
                return
            end
            P = obj.Pi;
            if isempty(obj.bgeos) || isempty(P)
                warning('Set bgeos and P0 fields.')
                return
            else
                roxy = obj.bgeos.density(P);
            end
            
            % Input parsing
            p = inputParser;
            p.FunctionName = mfilename;
            p.addParameter('axes', [], @(x)isscalar(x) && isgraphics(x, 'axes'))
            p.addParameter('plottype', 'line', @(x)isrow(x) && ischar(x))
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
            y = double(obj.rhoi - roxy);
            y = [y; y(end)];
            
            % Plot the lines (density in 1000 kg/m^3)
            if isequal(lower(pr.plottype), 'stairs')
                lh = stairs(x, y/1000);
            elseif isequal(lower(pr.plottype), 'line')
                lh = line(x, y/1000);
            else
                error('Unknown value of parameter plottype.') %#ok<ERTAG> 
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
                s = sprintf('$\\rho - \\rho_{\\mathrm{%s}}$ [1000 kg/m$^3$]',obj.bgeos.name);
                ylabel(s, 'fontsize', 12)
            else
                xlim('auto')
            end
            
            % Legend
            legend(ah, 'off')
            gh = legend(ah, 'show','location','ne');
            gh.FontSize = 11;
        end
        
        function [ah, lh, gh] = plot_Z_of_r(obj, varargin)
            % Plot Z(r) where r is mean radius.
            
            % Don't bother if uninitialized
            zvec = obj.zi;
            if isempty(zvec)
                warning('Level z value not found; set bgeos, fgeos, and P0 fields.')
                return
            end
            
            % Input parsing
            p = inputParser;
            p.FunctionName = mfilename;
            p.addParameter('axes', [], @(x)isscalar(x) && isgraphics(x, 'axes'))
            p.addParameter('plottype', 'line', @(x)isrow(x) && ischar(x))
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
            y = double([zvec; zvec(end)]);
            
            % Plot the lines
            if isequal(lower(pr.plottype), 'stairs')
                lh = stairs(x, y);
            elseif isequal(lower(pr.plottype), 'line')
                lh = line(x, y);
            else
                error('Unknown value of parameter plottype.') %#ok<ERTAG> 
            end
            lh.LineWidth = 2;
            if isempty(pr.axes)
                lh.Color = [0.31, 0.31, 0.31];
            end
            if isempty(obj.name)
                lh.DisplayName = sprintf('using %s',obj.fgeos.name);
            else
                lh.DisplayName = obj.name;
            end
            
            % Style and annotate axes
            if isempty(pr.axes)
                ah.Box = 'on';
                xlabel('Level surface radius, $s/s_0$', 'fontsize', 12)
                ylabel('$Z$', 'fontsize', 12)
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
                error('Unknown value of parameter plottype.') %#ok<ERTAG> 
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
        
        function [ah, lh] = plot_moi_contribution(obj, varargin)
            % Plot relative contribution to MoI by depth.
            
            p = inputParser;
            p.addParameter('axes',[],@(x)isscalar(x)&&isgraphics(x, 'axes'))
            p.addParameter('cumulative',false,@(x)isscalar(x)&&islogical(x))
            p.addParameter('noisecancel',false,@(x)isscalar(x)&&islogical(x))
            p.parse(varargin{:})
            pr = p.Results;
            
            % Prepare the data
            x = obj.si/obj.s0;
            if pr.cumulative
                y = obj.NMoI('csum');
                y = y/y(end);
            else
                y = obj.NMoI('none');
                if pr.noisecancel
                    y(x > 0.99) = nan;
                end
                y = y/max(abs(y));
            end
            
            % Prepare the canvas
            if isempty(pr.axes)
                fh = figure;
                set(fh, 'defaultTextInterpreter', 'latex')
                set(fh, 'defaultLegendInterpreter', 'latex')
                ah = axes;
            else
                ah = pr.axes;
                axes(ah)
            end
            hold(ah, 'on')
            
            % Plot the lines
            lh = plot(x, y, 'LineWidth',2);
            
            % Style and annotate axes
            if isempty(pr.axes)
                ah.Box = 'on';
                xlabel('Normalized level surface mean radius, $z_i=s_i/s_0$', 'fontsize', 12)
                if pr.cumulative
                    ylabel('$I(z>z_i)$ [normalized]', 'fontsize', 12)
                else
                    ylabel('$I(z_i)$ [normalized]', 'fontsize', 12)
                end
            end
        end
    end % End of visulaizers block
    
    %% Reporters/exporters
    methods (Access = public)
        function T = report_card(obj, obs)
            % REPORT_CARD Table summary of model's vital statistics.
            
            % Minimal checks
            narginchk(1,2);
            try
                obj.J2;
            catch
                warning('Uncooked object.') %#ok<*WNTAG>
                return
            end
            
            % Basic table
            vitals = {'Mass [kg]'; 'R_eq [km]'; 'J2'; 'J4'; 'J6'; 'J8'; 'NMoI'};
            TOF1 = [obj.M; obj.a0/1e3; obj.J2; obj.J4; obj.J6; obj.J8; obj.NMoI];
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
                oa0 = obs.a0/1e3;
            catch
                oa0 = NaN;
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
                oname = obs.name;
            catch
                oname = [];
            end
            OBS1 = [oM; oa0; oJ2; oJ4; oJ6; oJ8; oNMoI];
            OBS1 = double(OBS1);
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
            T = [T table(DIFF, 'VariableNames', {'frac_diff'})];
        end
        
        function s = to_struct(obj, rdc, keepss)
            % Convert object to static struct keeping only essential fields.
            
            if nargin < 2, rdc = 1; end % 0=none, 1=to double, 2=to single, 3=to scalars
            if nargin < 3, keepss = false; end % keep ss e.g. to help essample
            
            s.name   = obj.name;
            s.M      = obj.M;
            s.s0     = obj.s0;
            s.a0     = obj.a0;
            s.rhobar = obj.rhobar;
            s.mrot   = obj.mrot;
            s.qrot   = obj.qrot;
            s.J2     = obj.J2;
            s.J4     = obj.J4;
            s.J6     = obj.J6;
            s.J8     = obj.J8;
            if obj.opts.toforder == 7
                s.J10 = obj.Js(6);
                s.J12 = obj.Js(7);
                s.J14 = obj.Js(8);
            end
            s.NMoI   = obj.NMoI;
            s.si     = obj.si;
            s.ai     = obj.ai;
            s.rhoi   = obj.rhoi;
            s.Pi     = obj.Pi;
            s.mi     = obj.mi;
            
            if rdc > 0
                s = structfun(@double, s, 'UniformOutput', false);
                s.name = obj.name;
            end
            if rdc > 1
                s = structfun(@single, s, 'UniformOutput', false);
                s.name = obj.name;
            end
            if rdc > 2
                s.si     = [];
                s.ai     = [];
                s.rhoi   = [];
                s.Pi     = [];
                s.mi     = [];
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
            else
                s.ss = [];
            end
        end
        
        function T = to_table(obj)
            % Return a table of critical quantities.
            
            T = table;
            T.ai = double(obj.ai);
            T.si = double(obj.si);
            T.rhoi = double(obj.rhoi);
            T.Pi = double(obj.Pi);
            T.mi = double(obj.mi);
        end
        
        function to_ascii(obj, fname)
            % Export the state of the model as ascii file.
            
            % File name
            if nargin < 2
                fprintf('Usage:\n\ttof.to_ascii(filename)\n')
                return
            end
            validateattributes(fname, {'char'}, {'row'}, '', 'fname', 1)
            
            % Open file
            fid = fopen(fname,'wt');
            cleanup = onCleanup(@()fclose(fid));
            
            % Write the header
            fprintf(fid,'# Rotating fluid planet modeled by %dth-order Theory of Figures.\n',...
                        obj.opts.toforder);
            fprintf(fid,'#\n');
            fprintf(fid,'# Model name: %s\n', obj.name);
            fprintf(fid,'#\n');
            fprintf(fid,'# Scalar quantities:\n');
            fprintf(fid,'# N layers = %d\n',obj.N);
            fprintf(fid,'# Mass M = %g kg\n', obj.M);
            fprintf(fid,'# Mean radius       s0 = %0.6e m\n', obj.s0);
            fprintf(fid,'# Equatorial radius a0 = %0.6e m\n', obj.a0);
            fprintf(fid,'# Rotation period P = %0.6g s\n', obj.period);
            fprintf(fid,'# Rotation parameter m = %0.6f\n', obj.mrot);
            fprintf(fid,'# Rotation parameter q = %0.6f\n', obj.qrot);
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
            fprintf(fid,'# a_i   - equatorial radius of level surface i\n');
            fprintf(fid,'# rho_i - density on level surfaces i\n');
            fprintf(fid,'# P_i   - pressure on level surface i\n');
            fprintf(fid,'# m_i   - mass below level surface i\n');
            fprintf(fid,'#\n');
            
            % Write the data
            fprintf(fid,'# Column data:\n');
            fprintf(fid,'# %-4s  ','i');
            fprintf(fid,'%-10s  ','s_i','a_i');
            fprintf(fid,'%-7s  ','rho_i');
            fprintf(fid,'%-10s  ','P_i','m_i');
            fprintf(fid,'\n');
            for k=1:obj.N
                fprintf(fid,'  %-4d  ',k);
                fprintf(fid,'%10.4e  ', obj.si(k));
                fprintf(fid,'%10.4e  ', obj.ai(k));
                fprintf(fid,'%7.1f  ', obj.rhoi(k));
                fprintf(fid,'%10.4e  ', obj.Pi(k), obj.mi(k));
                fprintf(fid,'\n');
            end
        end
    end % End of reporters/exporters block
    
    %% Private (or obsolete) methods
    methods (Access = private)
        function val = Ki(obj)
            P = obj.Pi;
            ro = obj.rhoi;
            if isempty(P) || isempty(ro)
                val = [];
            else
                dPdro = sdderiv(ro,P);
                val = ro.*dPdro;
            end
        end
        
        function y = P_mid(obj)
            % Pressure interpolated to halfway between level surfaces.
            
            v = double(obj.Pi);
            if isempty(v), y = []; return, end
            x = double(obj.si);
            xq = [(x(1:end-1) + x(2:end))/2; x(end)/2];
            y = interp1(x, v, xq, 'pchip');
        end
        
        function val = mzi(obj)
            % heavy element mass below level i
            z = obj.zi;
            if isempty(obj.si) || isempty(obj.rhoi) || isempty(z)
                val = [];
            else
                rho = obj.rhoi;
                s = obj.si;
                val(obj.N) = 4*pi/3*rho(obj.N)*s(obj.N)^3*z(obj.N);
                for k=obj.N-1:-1:1
                    cz = max(min(z(k), 1), 0);
                    val(k) = val(k+1) + 4*pi/3*rho(k)*(s(k)^3 - s(k+1)^3)*cz;
                end
                val = val';
            end
        end
        
        function val = zi(obj)
            % heavy element mass fraction on level i
            P = obj.Pi;
            if isempty(obj.bgeos) || isempty(obj.fgeos) || isempty(P)
                val = [];
            else
                roxy = obj.bgeos.density(P);
                roz = obj.fgeos.density(P);
                ro = obj.rhoi;
                val = (1./ro - 1./roxy)./(1./roz - 1./roxy);
                val(~isfinite(val)) = 0;
            end
        end
        
        function val = M_Z(obj)
            try
                val = obj.mzi(1);
            catch
                val = [];
            end
        end
        
    end % End of private methods block
    
    %% Access and pseudo-access methods
    methods
        function set.name(obj,val)
            if ~isempty(val)
                validateattributes(val, {'char'}, {'row'})
            end
            obj.name = val;
        end
        
        function set.mass(obj,val)
            validateattributes(val,{'numeric'},{'positive','scalar'})
            obj.mass = val;
        end
        
        function set.radius(obj,val)
            validateattributes(val,{'numeric'},{'positive','scalar'})
            obj.radius = val;
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
        
        function set.eos(obj,val)
            if isempty(val)
                obj.eos = [];
                return
            end
            if ~isa(val,'barotropes.Barotrope')
                error('eos must be a valid instance of class Barotrope') %#ok<ERTAG> 
            end
            obj.eos = val(:);
        end
        
        function set.bgeos(obj,val)
            if isempty(val)
                obj.bgeos = [];
                return
            end
            if ~isa(val,'barotropes.Barotrope')
                error('bgeos must be a valid instance of class Barotrope') %#ok<ERTAG> 
            end
            obj.bgeos = val;
        end
        
        function set.fgeos(obj,val)
            if isempty(val)
                obj.fgeos = [];
                return
            end
            if ~isa(val,'barotropes.Barotrope')
                error('fgeos must be a valid instance of class Barotrope') %#ok<ERTAG> 
            end
            obj.fgeos = val;
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
                    length(obj.si),length(obj.rhoi)) %#ok<ERTAG> 
            end
        end
        
        function val = Ui(obj, ind)
            % Following Nettelmann 2017 eqs. B3 and B.4, assuming equipotential.
            if isempty(obj.A0), val = []; return, end
            
            val = -obj.G*obj.mass/obj.s0^3*obj.si.^2.*obj.A0;
            if nargin > 1
                val = val(ind);
            end
        end
        
        function val = Pi(obj, ind)
            if isempty(obj.Ui) || isempty(obj.rhoi) || isempty(obj.P0)
                val = [];
                return
            end
            n = obj.N;
            U = obj.Ui;
            rho = obj.rhoi;
            r = obj.si;
            gradU = zeros(n,1);
            val = zeros(n,1);
            gradU(1) = (U(1) - U(2))/(r(1) - r(2));
            gradU(2:n-1) = (U(1:n-2) - U(3:n))./(r(1:n-2) - r(3:n));
            gradU(n) = (U(n-1) - U(n))/(r(n-1) - r(n));
            intgrnd = rho.*gradU;
            val(1) = obj.P0;
            switch lower(obj.opts.prsmeth)
                case 'eq6'
                    for k=1:n-1
                        val(k+1) = val(k) + 0.5*(rho(k) + rho(k+1))*(U(k+1) - U(k));
                    end
                case 'trapz'
                    for k=1:n-1
                        val(k+1) = val(k) + 0.5*(r(k) - r(k+1))*(intgrnd(k) + intgrnd(k+1));
                    end
                case 'simps'
                    h = mean(abs(diff(r)));
                    val(2) = val(1) + (h/2)*(intgrnd(1) + intgrnd(2)); % trapz prime
                    for k=3:n
                        val(k) = val(k-2) + ...
                            (h/3)*(intgrnd(k-2) + 4*intgrnd(k-1) + intgrnd(k));
                    end
                otherwise
                    error('Unknown pressure integral method.') %#ok<ERTAG> 
            end
            if nargin > 1
                val = val(ind);
            end
        end
        
        function val = get.M(obj)
            if isempty(obj.si) || isempty(obj.rhoi)
                val = [];
            else
                val = obj.total_mass(obj.opts.masmeth);
            end
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
                end % TODO: switch to cumtrapz
                val = val';
            end
        end
        
        function val = get.rhobar(obj)
            if isempty(obj.M) || isempty(obj.si)
                val = [];
            else
                val = obj.M/(4*pi/3*obj.s0^3);
            end
        end
        
        function val = get.wrot(obj)
            val = 2*pi./obj.period;
        end
        
        function val = get.qrot(obj)
            GM = obj.G*obj.mass;
            val = obj.wrot^2.*obj.radius^3./GM;
        end
        
        function val = get.mrot(obj)
            GM = obj.G*obj.mass;
            val = obj.wrot^2.*obj.s0^3./GM;
        end
        
        function set.mrot(~,~)
            error('TOFPLANET:deprecation',...
                'Setting mrot is deprecated; set a reference period instead.')
        end
        
        function val = get.a0(obj)
            if isempty(obj.si)
                val = [];
            else
                val = obj.aos*obj.s0;
            end
        end
        
        function val = get.ai(obj)
            if isempty(obj.si) || isempty(obj.ss)
                val = [];
            else
                val = ones(size(obj.si));
                for k=1:obj.N
                    shp = obj.ss.s0(k)*Pn(0,0) + obj.ss.s2(k)*Pn(2,0) + ...
                        obj.ss.s4(k)*Pn(4,0) + obj.ss.s6(k)*Pn(6,0) + ...
                        obj.ss.s8(k)*Pn(8,0);
                    val(k) = obj.si(k)*(1 + shp);
                end
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
        
        function val = get.J10(obj)
            if isempty(obj.Js) || (length(obj.Js) < 6)
                val = 0;
            else
                val = obj.Js(6);
            end
        end
        
        function val = get.J12(obj)
            if isempty(obj.Js) || (length(obj.Js) < 7)
                val = 0;
            else
                val = obj.Js(7);
            end
        end
        
    end % End of access methods block
    
    %% Static methods
    methods (Static)
        function r = kth_surface(k,ss)
            % Return @r(theta) for level surface k.
            
            shp = @(mu) ss.s0(k)*Pn(0,mu) + ss.s2(k)*Pn(2,mu) + ...
                ss.s4(k)*Pn(4,mu) + ss.s6(k)*Pn(6,mu) + ...
                ss.s8(k)*Pn(8,mu);
            r = @(teta) 1 + shp(cos(teta));
        end
    end % End of static methods block
end % End of classdef

%% Class-related functions
function [x,w] = gauleg(x1,x2,n)
%GAULEG Calculate abscissas and weights for Gauss-Legendre n-point quadrature.
%   [x,w] = GAULEG(x1,x2,n) returns the abscissas x and weights w that can be
%   used to evaluate the definite integral, I, of a function well approximated
%   by an (2n - 1) degree polynomial in the interval [x1,x2] using the
%   Gauss-Legendre formula:
%
%       I = sum(w.*f(x))
%
%   Algorithm
%     This function is based on the C++ implementation of a routine with the
%     same name in Numerical Recipes, 3rd Edition. But in several places I opt
%     for readability over performance, on the assumption that this function is
%     most likely to be called in a setup routine rather than in an inner-loop
%     computation.
%
%   Example
%     fun = @(x)sin(x);
%     [x,w] = gauleg(0,pi,6);
%     I_adaptive = integral(fun,0,pi)
%     I_gaussleg = sum(w.*fun(x))
%
% Author: Naor Movshovitz (nmovshov at google dot com)
%         Earth and Planetary Sciences, UC Santa Cruz
%
% Reference: William H. Press, Saul A. Teukolsky, William T. Vetterling, and
% Brian P. Flannery. 2007. Numerical Recipes 3rd Edition: The Art of Scientific
% Computing (3 ed.). Cambridge University Press, New York, NY, USA.

% Input parsing and minimal assertions
narginchk(3,3)
nargoutchk(2,2)
validateattributes(x1,{'numeric'},{'scalar','finite','real'},1)
validateattributes(x2,{'numeric'},{'scalar','finite','real'},2)
validateattributes(n,{'numeric'},{'scalar','finite','integer','>=',2},3)
assert(x2 > x1, 'Interval must be positive.');

% Local variables
tol = 1e-14;
m = ceil(n/2);
xmid = (x1 + x2)/2;
dx = (x2 - x1);
x = NaN(1,n);
w = NaN(1,n);

% Main loop
for j=1:m
    % Get j-th root of Legendre polynomial Pn, along with Pn' value there.
    z = cos(pi*((j - 1) + 0.75)/(n + 0.5)); % initial guess for j-th root
    while true
        % Calculate Pn(z) and Pn-1(z) and Pn'(z)
        p = NaN(1,n+1);
        p(1) = 1;
        p(2) = z;
        for k=2:n
            pkm1 = p(k);
            pkm2 = p(k-1);
            pk = (1/k)*((2*k - 1)*z*pkm1 - (k - 1)*pkm2);
            p(k+1) = pk;
        end
        pn = p(end);
        pp = (n*p(end-1) - n*z*p(end))/(1 - z^2);
        
        % And now Newton's method (we are hopefully very near j-th root)
        oldz = z;
        z = z - pn/pp;
        if abs(z - oldz) < tol, break, end
    end
    
    % Now use j-th root to get 2 abscissas and weights
    x(j)     = xmid - z*dx/2; % Scaled abscissa left of center
    x(n+1-j) = xmid + z*dx/2; % Scaled abscissa right of center
    w(j)     = dx/((1 - z^2)*pp^2);
    w(n+1-j) = w(j);
end

% Verify and return
assert(all(isfinite(x)))
assert(all(isfinite(w)))
end

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
