        function [rcore, mcore, icore] = locate_core(obj, how)
            % Look for a solid core.
            %
            % [rcore, mcore, icore] = LOCATE_CORE(how) returns the estimated
            % radius, mass, and starting level index of the planet's core. There
            % are several options for how to decide what part of the planet, if
            % any, is a core, and these are specified using the string argument
            % 'how', in the following ways:
            %
            % LOCATE_CORE('byeos') locates the innermost contiguous block of
            % levels associated with the same Barotrope object. This gives the
            % correct answer, by definition, in the rare cases where the planet
            % was constructed with a particular barotrope-based model in mind, and
            % relaxed using obj.relax_to_barotrope instead of obj.relax_to_HE.
            %
            % LOCATE_CORE('byrho') attempts to identify the innermost block of
            % contiguous levels following a detectable density jump. This block is
            % considered to be a core if its density is higher than obj.rhoc.
            %
            % LOCATE_CORE('byZ') attempts to identify the innermost block of
            % contiguous levels following a detectable jump in heavy elements
            % abundance (obj.zi). This jump is considered to define a core if it
            % leads to z greater than obj.zeec.
            %
            % LOCATE_CORE('byzeec') returns the mass below the outermost layer in
            % the inner half of the planet where obj.zi > obj.zeec.
            %
            % If a core cannot be detected for any reason the function returns
            % [0,0,0].
            
            if nargin == 1 && nargout == 0
                help TOFPlanet.locate_core
                return
            end
            how = validatestring(lower(how),{'byeos','byrho','byz','byzeec'});
            
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
                        mcore = 0;
                        icore = 0;
                    else
                        if zlay, cind = cind + 1; end
                        icore = cind;
                        rcore = obj.si(icore);
                        mcore = obj.mi(icore);
                    end
                case 'byrho'
                    % Define core using a density jump
                    dvec = double(obj.rhoi/obj.rhobar)';
                    deltas = diff(dvec);
                    cind = peakfinder(deltas);
                    if isempty(cind)
                        rcore = 0;
                        mcore = 0;
                        icore = 0;
                    else
                        icore = cind(end)+1;
                        rcore = obj.si(icore);
                        mcore = obj.mi(icore);
                        if isempty(obj.rhoc)
                            error('obj.rhoc not set; set to zero to ignore.')
                        end
                        if obj.rhoi(icore) < obj.rhoc
                            icore = 0;
                            rcore = 0;
                            mcore = 0;
                        end
                    end
                case 'byzeec'
                    % Define a dilute core using r/R < 0.5 & zi > zeec
                    rvec = double(obj.si/obj.s0);
                    zvec = double(obj.zi);
                    if isempty(obj.zeec)
                        error('obj.zeec not set; set to zero to ignore.')
                    end
                    cind = find((rvec < 0.5) & (zvec > obj.zeec), 1, 'first');
                    if isempty(cind)
                        icore = 0;
                        rcore = 0;
                        mcore = 0;
                    else
                        icore = cind;
                        rcore = obj.si(icore);
                        mcore = obj.si(icore);
                    end
                case 'byz'
                    % Define core using Z jump
                    zvec = double(obj.zi);
                    zvec(~isfinite(zvec)) = 0;
                    deltas = diff(zvec);
                    cind = peakfinder(deltas(2:end)); % first layer often noisy
                    if isempty(cind)
                        icore = 0;
                        rcore = 0;
                        mcore = 0;
                    else
                        icore = cind(end)+2; % remember deltas "lost" 2 layers
                        rcore = obj.si(icore);
                        mcore = obj.mi(icore);
                        if isempty(obj.zeec)
                            error('obj.zeec not set; set to zero to ignore.')
                        end
                        if obj.zi(icore) < obj.zeec
                            icore = 0;
                            rcore = 0;
                            mcore = 0;
                        end
                    end
            end
        end
