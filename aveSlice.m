classdef aveSlice < kCut
    % MEANSLICE Contains the 2D (spanwise averaged) mean flow

        
    properties
        p0in;
        Uinf;
        roinf;
        muinf;
        Tinf;
        T0in;
        pin;
        dsdyThresh;% = -50;
        use_unsflo = false;
        blEdgeMode;
        blEdgeStore = [];
        iSmoothBL =  false;
        nSmooth;
        label;
    end

    properties (Dependent = true, Hidden = true)
        Msurf;          % Surface Mach No
        Psurf;
        theta;          % Momentum thickness
        theta_k;
        thetaStar;      % K.E. Thickness
        H;              % Shape factor
        H_k;            % Kinematic shape factor
        H_k_Whitfield;  % Whitfield correlation for Hk for adiabatic flows
        H_ke;           % K.E. shape factor = thetaStar/theta
        H_rho;          % Density shape factor H** = delta**/theta
        H1;             
        Us;             % Slip velocity
        %delta99;       % BL thickness
        %delta99_unsflo;
        delta;          % Approx BL thickness
        delStar;        % Displacemnt thickness
        delta99_k;      % BL thickness, incompressible definition
        delStar_k;      % Displacemnt thickness, incompressible definition
        delRho;         % Density thickness, delta**
        dsdy;           % Wall normal entropy gradient
        %BLedgeInd;      % j index of detected BL edge
        U;              % Wall-parallel velocity
        Ue;             % BL edge velocity
        Me;             % BL edge Mach number
        Res;            % Surface distance Reynolds No
        Pr_nondim;
        blPr;           % Componant of cd due to production of tke
        blPr_dimensional;
        blPr_eq;        % Coles' equilibrium prodution
        cPr;            % Non-dimensional local production
        cd_inner;       % Componant of cd dow to mean strain
        cd_mean_strain;
        cd;             % Dissipation coefficient
        tau_w;          % Wall shear stress
        u_tau;          % Friction velocity
        dUdy;
        dTdy;
        cf;
        ct;             % ctau based on Pr/(1-Us);
        ctau;
        ctau_max;
        Re_theta;
        iPS;
        iMaxPr;
        iEq;
        Re_theta_ps;
        pdyn            % Dynamic pressure, 0.5*rho*V^2
        nu_e            % Boundary layer edge viscosity
        alpha1;
        alpha2;
        zeta;           % Entropy loss coefficient
        Yp;             % Stagnation pressure loss coefficient
        p0out;
        yplus;
        wall_scale;     % Coordinade for wall distance = 1;
        mut;            % Turbulent viscosity
        mut_ratio;
        RT;             % Hydra stupe turbulent Re
%         wallDist;
        pr_aligned_s;
        pr_aligned_n;
        pr_aligned_t;
    end

    methods
        function obj = aveSlice(blk, gas, bcs, iPS)
            if nargin < 4
                iPS = false;
            end
            obj@kCut(blk, gas, bcs, iPS);
            disp('Constructing aveSlice')
            obj.nSmooth = round(obj.niBL/50);
        end

        function getBCs(obj, inlet_blocks, is)
            if nargin < 3
                is = 40:100;
            end
            if nargin < 2 || isempty(inlet_blocks)
                inlet_blocks = obj.blk.inlet_blocks{1};
            end
            Mnow = obj.M;
            Unow = obj.vel;
            ronow = obj.ro;
            munow = obj.mu;
            %Mnow = Mnow{inlet_blocks};
            pnow = obj.p;
            Tnow = obj.T;
            %pnow = pnow{inlet_blocks};
            
            p0 = [];
            Uinf = [];
            muinf = [];
            roinf = [];
            Tinf = [];
            pinf = [];
            for i=1:length(inlet_blocks)
                nj = obj.blk.blockdims(inlet_blocks(i),2);
                js = 1:nj;
                if obj.blk.next_patch{inlet_blocks(i)}.jm == 3 && ...
                        obj.blk.next_block{inlet_blocks(i)}.jm == 0
                    inds = obj.BLedgeInd;
                    ind = max(inds(is));
                    js = floor(1.2*ind):nj;
                elseif obj.blk.next_patch{inlet_blocks(i)}.jp == 3 && ...
                        obj.blk.next_block{inlet_blocks(i)}.jp == 0
                    inds = obj.BLedgeInd;
                    ind = max(inds(is));
                    js = (nj-floor(1.2*ind)):nj;
                end
                p0now = pnow{inlet_blocks(i)}.*(1+((obj.gas.gam - 1)/2)*Mnow{inlet_blocks(i)}.^2).^(obj.gas.gam/(obj.gas.gam-1));
                p0 = [p0 p0now(is,js)];
                Tinf = [Tinf Tnow{inlet_blocks(i)}(is,js)];
                Uinf = [Uinf Unow{inlet_blocks(i)}(is,js)];
                muinf = [muinf munow{inlet_blocks(i)}(is,js)];
                roinf = [roinf ronow{inlet_blocks(i)}(is,js)];
                pinf = [pinf pnow{inlet_blocks(i)}(is,js)];
            end
            obj.p0in = mean(p0,'all');
            obj.Tinf = mean(Tinf, 'all');
            obj.Uinf = mean(Uinf,'all');
            obj.pin = mean(pinf, 'all');
            obj.muinf = mean(muinf,'all');
            obj.roinf = mean(roinf,'all');
            obj.T0in = obj.Tinf+obj.Uinf^2/(2*obj.gas.cp);
%             obj.wallDist = obj.blk.y;
        end

        function value = smooth_bl_edge(obj, x, y, xSafe)

            if nargin < 4
                xSafe = 0.25;
            end
            ilast = obj.x2ind(xSafe);
            points2sample = zeros(size(x));

            for i = ilast:length(x)
                if abs((y(i)-y(ilast))/y(ilast)) > 0.2
                    points2sample(i) = true;
                else
                    ilast = i;
                end
            end

            %points2sample = abs(diff(y)./y(1:end-1))>0.1;
            ytmp = y(~points2sample);
            xtmp = x(~points2sample);
            for i = find(points2sample)
                y(i) = round(interp1(xtmp,ytmp,x(i),'linear'));
            end
            value = y;
        end

        function value = smooth_dist(obj, prop, nsmooth)

            if nargin < 3 || isempty(nsmooth)
                nsmooth = obj.nSmooth;
            end
            
            nans = isnan(prop);
            value = smooth(prop, nsmooth, 'rloess');
            value(nans) = NaN;

        end

        function value = get.dsdy(obj)
            %s = obj.oGridProp('s');

            %value = (s(:,2:end)-s(:,1:end-1))./(obj.yBL(:,2:end)-obj.yBL(:,1:end-1));
            value = obj.blNormGrad('s');
        end

        function value = get.dUdy(obj)
%             u = obj.oGridProp('U');
%             value = (u(:,2:end)-u(:,1:end-1))./(obj.yBL(:,2:end)-obj.yBL(:,1:end-1));
            value = obj.blNormGrad('U');
        end 

        function value = get.dTdy(obj)
            t = obj.oGridProp('T');
            value = (t(:,2:end)-t(:,1:end-1))./(obj.yBL(:,2:end)-obj.yBL(:,1:end-1));
            value = obj.blNormGrad('T');
        end 

        function value = get.Msurf(obj)
            disp('Calculating surface M')
            pnow = obj.oGridProp('p');
            value = sqrt((2/(obj.gas.gam - 1)) * ( (pnow(:,1)/obj.p0in).^(-(obj.gas.gam-1)/obj.gas.gam) - 1));
        end

        function value = get.Psurf(obj)
            disp('Calculating surface p')
            pnow = obj.oGridProp('p');
            value = pnow(:,1);
        end
       

        function value = BLedgeInd(obj, mode, debug)
            if nargin < 2 || isempty(mode)
                if ~isempty(obj.blEdgeStore)
                    value = obj.blEdgeStore;
                    return
                else
                    if ~isempty(obj.blEdgeMode)
                        mode = obj.blEdgeMode;
                    else
                        mode = "sCombined";
                    end
                end
            end
            if nargin < 3 || isempty(debug)
                debug = false;
            end
            
            fprintf('BL edge detection mode: %s\n',mode)
            switch mode
                case "sGradThresh"
                    temp = obj.dsdy;
                    for i=1:size(temp,1)
                        j=1;
                        if isempty(obj.dsdyThresh)
                            thresh = 0.02*temp(i,1);% + 0.98*(temp(i,ceil(size(temp,2)/2))-temp(i,1));
                        else
                            thresh = obj.dsdyThresh;
                        end
                        while temp(i,j) < thresh && j<size(temp,2)
                            j = j+1;
                        end
                        value(i) = min(j+1,size(temp,2));
                        [~, value(i)] = max(u(i,1:value(i)));
                    end
                case "sThresh"
                    temp = obj.oGridProp('s');
                    for i=1:size(temp,1)
                        if i==1255
                            disp('')
                        end
%                         j=size(temp,2);
                        j=1;
                        se = temp(i,end);
                        sw = temp(i,1);
                        se;
                        while temp(i,j) > se + 0.02*(sw-se) && j<size(temp,2)
                            j = j+1;
                        end
                        value(i) = min(j+1,size(temp,2));
                        %[~, value(i)] = max(u(i,1:value(i)));
                    end
                case "sCombined"
                    snow = obj.oGridProp('s');
                    dsdy = obj.dsdy;
                    for i=1:size(snow,1)

                        if i == 4094
                            disp('');
                        end
                        
                        if isempty(obj.dsdyThresh)
                            dsdynow = smooth(dsdy(i,:), round(size(dsdy,2)/20));
                            % dsdythresh = -0.01*max(abs(dsdy(i,:)));% + 0.98*(temp(i,ceil(size(temp,2)/2))-temp(i,1));
                            dsdythresh = -0.01*max(abs(dsdynow));
                        else
                            dsdythresh = obj.dsdyThresh;
                        end

                        inds = max(find(abs(dsdy(i,:)) > abs(dsdythresh)));
                        inds = max(inds):size(dsdy,2);
                        % inds = abs(dsdy(i, :)) < abs(dsdythresh);
                        splateau = mean(snow(i,inds));
                        sthresh  = 0.98*splateau + 0.02*max(abs(snow(i,:)));

                        j=1;
                        while ((snow(i,j) > sthresh) ...% && snow(i,j) > snow(i, j+1)) ...
                                || dsdy(i,j) < dsdythresh) ...
                                && j<size(snow,2)
                            j = j+1;
                        end
                        value(i) = min(j+1,size(snow,2));
%                         if value(i) == 2
%                             value(i) = NaN;
%                         end
                        
                    end
                case "unsflo"
                    del99 = obj.delta99("unsflo");
                    for i=1:length(del99)
                        ynow = obj.yBL(i,:);
                        [~, value(i)] = min(abs(del99(i) - ynow));
                        [~, value(i)] = max(u(i,1:value(i)));
                        %value(i) = min(size(temp,2), indnow+1);
                    end
                case "sInteg"
                    u = obj.oGridProp('U');
                    del99 = obj.delta99("sInteg");
                    for i=1:length(del99)
                        ynow = obj.yBL(i,:);
                        [~, value(i)] = min(abs(del99(i) - ynow));
                        [~, value(i)] = max(u(i,1:value(i)));
                        %value(i) = min(size(temp,2), indnow+1);
                    end
                case "p0Integ"
                    del99 = obj.delta99("p0Integ");
                    for i=1:length(del99)
                        ynow = obj.yBL(i,:);
                        [~, value(i)] = min(abs(del99(i) - ynow));
                        %value(i) = min(size(temp,2), indnow+1);
                    end

                case "integralCorrelation"
                    mode = obj.blEdgeMode;
                    obj.blEdgeMode = "sCombined";
                    
                    U = obj.U;
                    del = obj.delta99;
                    inds = obj.thickness2ind(del);
                    for i=1:length(del)
                        Uenow = 0.99*U(i, inds(i));
                        [~, value(i)] = min(abs(U(i,1:inds(i)) - Uenow));
                    end
                    obj.blEdgeMode = mode;

            end
            value(value<5) = nan;
            if obj.iSmoothBL
                value = round(obj.smooth_dist(value));
            end

            obj.blEdgeStore = value;
            
        end

        function [xedge, yedge] = getBLedgeCoords(obj,mode)
            if nargin < 2 || isempty(mode)
                if ~isempty(obj.blEdgeMode)
                    mode = obj.blEdgeMode;
                else
                    mode = "sCombined";
                end
            end
            inds = obj.BLedgeInd(mode);
            for i=1:length(inds)
                xedge(i) = obj.xO(i,inds(i));
                yedge(i) = obj.yO(i,inds(i));
            end
        end

        function value = get.iPS(obj)
            x = obj.xSurf;
            M = obj.smooth_dist(obj.Msurf);
            pr = obj.smooth_dist(obj.blPr);

            i = length(x);
            while M(i) < 1.02
                i = i-1;
            end
            prLast = inf;
            while pr(i) < prLast
                prLast = pr(i);
                i = i-1;
            end
            value = i;
        end

        function value = get.iMaxPr(obj)
            [~, value] = max(obj.blPr);
        end

        function value = get.iEq(obj)
            i = 20;
            pr = obj.smooth_dist(obj.blPr);
            while pr(i) > pr(i-1)
                i = i+1;
            end
            value = i;
        end

        function value = delta99_unsflo(obj)
            %inds = obj.BLedgeInd;
            
        end

        function value = delta99(obj, mode)
            if nargin < 2 || isempty(mode)
                if ~isempty(obj.blEdgeMode)
                    mode = obj.blEdgeMode;
                else
                    mode = "sCombined";
                end
            end
            switch mode
                case "sGradThresh"
                    inds = obj.BLedgeInd("sGradThresh");
                    value(isnan(inds)) = nan;
                    for i=find(~isnan(inds))
                        value(i) = obj.yBL(i,inds(i));
                    end
                case "unsflo"
                    ugrad = abs(obj.dUdy);
                    tgrad = abs(obj.dTdy);
                    for i=1:size(obj.yBL,1)
                        ys = obj.yBL(i,:);
                        %value(i) = obj.yBL(i,inds(i));
                        int1 = trapz(ys, ys.*(ugrad(i,:)));%+tgrad(i,:)));
                        int2 = trapz(ys, (ugrad(i,:)));%+tgrad(i,:)));
                        value(i) = 2.5*int1/int2;
                    end
                case "sInteg"
                    sgrad = abs(obj.blNormGrad('s'));
                    for i=1:size(obj.yBL,1)
                        ys = obj.yBL(i,:);
                        %value(i) = obj.yBL(i,inds(i));
                        int1 = trapz(ys, ys.*(sgrad(i,:)));%+tgrad(i,:)));
                        int2 = trapz(ys, (sgrad(i,:)));%+tgrad(i,:)));
                        value(i) = 2.5*int1/int2;
                    end
                case "p0Integ"
                    p0grad = abs(obj.blNormGrad('p0'));
                    for i=1:size(obj.yBL,1)
                        ys = obj.yBL(i,:);
                        %value(i) = obj.yBL(i,inds(i));
                        int1 = trapz(ys, ys.*(p0grad(i,:)));%+tgrad(i,:)));
                        int2 = trapz(ys, (p0grad(i,:)));%+tgrad(i,:)));
                        value(i) = 2.5*int1/int2;
                    end
                case "sThresh"
                    inds = obj.BLedgeInd("sThresh");
                    value(isnan(inds)) = nan;
                    for i=find(~isnan(inds))
                        value(i) = obj.yBL(i,inds(i));
                    end
                case "sCombined"
                    inds = obj.BLedgeInd("sCombined");

                    value(isnan(inds)) = nan;
                    for i=find(~isnan(inds))
                        value(i) = obj.yBL(i,inds(i));
                    end
                case "integralCorrelation"
                    inds = obj.BLedgeInd("integralCorrelation");
                    value(isnan(inds)) = nan;
                    for i=find(~isnan(inds))
                        value(i) = obj.yBL(i,inds(i));
                    end
            end
            value = obj.smooth_dist(value);
            
        end


        function value = get.U(obj)
            disp('Calculating U')
            unow = obj.oGridProp('u');
            vnow = obj.oGridProp('v');
            nnow = obj.n;
            value = zeros(size(obj.yBL));
            for i=1:size(obj.yBL,1)
                tang = [nnow(2,i); -nnow(1,i)];
                for j=1:size(obj.yBL,2)
                    velnow = [unow(i,j); vnow(i,j)];
                    % value(i,j) = -dot(tang, velnow - nnow(:,i)*dot(nnow(:,i),velnow));
                    value(i,j) = dot(velnow - dot(velnow,nnow(:,i)), tang);
                end
            end
        end

        function value = get.delStar(obj)
            inds = obj.BLedgeInd;
            ronow = obj.oGridProp('ro');
            Unow = obj.U;
            value = zeros(1,length(inds));
            value(isnan(inds)) = NaN;
            for i=find(~isnan(inds)) 
                roprof = ronow(i,1:inds(i));
                Uprof = Unow(i,1:inds(i));
                ro0 = ronow(i,inds(i));
                U0 = Unow(i,inds(i));
                integrand = 1 - ((roprof.*Uprof)./(ro0*U0));
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, integrand);
            end
            value = obj.smooth_dist(value);
        end

        function value = get.delStar_k(obj)
            inds = obj.BLedgeInd;
            Unow = obj.U;
            value = zeros(1,length(inds));
            value(isnan(inds)) = NaN;
            for i=find(~isnan(inds)) 
                integrand = 1 - Unow(i,1:inds(i))./Unow(i,inds(i));
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, integrand);
            end
            value = reshape(value, [], 1);
            value = obj.smooth_dist(value);
        end

        function value = get.delRho(obj)
            inds = obj.BLedgeInd;
            ronow = obj.oGridProp('ro');
            Unow = obj.U;
            value(isnan(inds)) = NaN;
            for i=find(~isnan(inds)) 
                roprof = ronow(i,1:inds(i));
                Uprof = Unow(i,1:inds(i));
                ro0 = ronow(i,inds(i));
                U0 = Unow(i,inds(i));
                integrand = (1 - roprof/ro0).*Uprof/U0;
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, integrand);
            end
            value = obj.smooth_dist(value);
        end

        function value = get.theta(obj)
            inds = obj.BLedgeInd;
            ronow = obj.oGridProp('ro');
            Unow = obj.U;
            value(isnan(inds)) = NaN;
            for i=find(~isnan(inds)) 
                roprof = ronow(i,1:inds(i));
                Uprof = Unow(i,1:inds(i));
                ro0 = ronow(i,inds(i));
                U0 = Unow(i,inds(i));
                integrand = (roprof.*Uprof/(ro0*U0)).*(1-Uprof/U0);
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, integrand);
            end
            value = obj.smooth_dist(value);
        end

        function value = get.theta_k(obj)
            inds = obj.BLedgeInd;
            Unow = obj.U;
            value(isnan(inds)) = NaN;
            for i=find(~isnan(inds)) 
                Uprof = Unow(i,1:inds(i));
                U0 = Unow(i,inds(i));
                integrand = (Uprof/U0).*(1-Uprof/U0);
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, integrand);
            end
            value = reshape(value, [], 1);
            value = obj.smooth_dist(value);
        end

        function value = get.thetaStar(obj)
            inds = obj.BLedgeInd;
            ronow = obj.oGridProp('ro');
            Unow = obj.U;
            value(isnan(inds)) = NaN;
            for i=find(~isnan(inds)) 
                roprof = ronow(i,1:inds(i));
                Uprof = Unow(i,1:inds(i));
                ro0 = ronow(i,inds(i));
                U0 = Unow(i,inds(i));
                integrand = (roprof.*Uprof/(ro0*U0)).*(1-Uprof.^2/U0^2);
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, integrand);
            end
            value = obj.smooth_dist(value);
        end

        function value = get.H(obj)
            value = obj.delStar./obj.theta;
        end

        function value = get.H_k(obj)
            value = obj.delStar_k./obj.theta_k;
        end

        function value = get.H_k_Whitfield(obj)
            H = obj.H;
            Me = obj.Msurf;
            value = (H-0.290*Me.^2)./(1+0.113*Me.^2);
        end

        function value = get.H_ke(obj)
            value = obj.thetaStar./obj.theta;
        end

        function value = get.H_rho(obj)
           value = obj.delRho./obj.theta;
        end

        function value = get.H1(obj)
            inds = obj.BLedgeInd;
            ronow = obj.oGridProp('ro');
            Unow = obj.U;
            value(isnan(inds)) = NaN;
            for i=find(~isnan(inds)) 
                roprof = ronow(i,1:inds(i));
                Uprof = Unow(i,1:inds(i));
                ro0 = ronow(i,inds(i));
                U0 = Unow(i,inds(i));
                integrand = (roprof.*Uprof)/(ro0*U0);
                ys = obj.yBL(i,1:inds(i));
                num(i) = trapz(ys, integrand);
            end
            value = num(i)./obj.theta;
        end

        function value = get.Us(obj)
            value = MISES_correlations.fUs(obj.H_ke, obj.H_k, obj.H);
        end

        function value = get.delta(obj)
            value = (obj.H + obj.H1).*obj.theta;
        end

        function value = get.Res(obj)
            value = obj.ssurf*obj.Uinf*obj.roinf/obj.muinf;
        end

        function value = obj.dUedx_eq(obj)

            A = 6.7;
            B = 0.75;

            Hnow = obj.H;
            Hknow = obj.Hk;
            thetanow = obj.theta;

            value = (1/(B*Hnow.*thetanow))*(0.5*obj.cf - ((Hknow-1)/(A*Hk))^2);
        end

        function value = ctau_eq(obj)

            A = 6.7;
            B = 0.75;
            Hknow = obj.H_k;
            
            value = (0.5/(A*A*B))*obj.H_ke.*(Hknow-1).^3./(obj.H.*Hknow.^2);

        end

        function value = cd_eq(obj)
            % 
            % Us = obj.Us;
            % value = 0.5*obj.cf'.*Us + obj.ctau_eq.*(1-Us);
            % 

            value = obj.cd_mean_strain + obj.blPr_eq;
        end


        function plt = blDevPlot(obj, prop, varargin) % ax, lims, xrange, fmt)

            x0 = min(obj.xSurf);
            x1 = max(obj.xSurf);

            defaultAx = gca;
            defaultLims = 'auto';
            defaultLineWidth = 1.5;
            defaultFontSize = 12;
            defaultFmtString = '';
            defaultXRange = [x0-1 x1+1];

            p = inputParser;

%             addRequired(p, 'prop');
            addParameter(p, 'ax', defaultAx);
            addParameter(p, 'lims', defaultLims);
            addParameter(p, 'xrange', defaultXRange);
            addParameter(p, 'fmt', defaultFmtString);
            addParameter(p, 'LineWidth', defaultLineWidth);
            addParameter(p, 'FontSize', defaultFontSize)

            parse(p, varargin{:})

%             if nargin < 3 || isempty(ax)
%                 ax = gca;
%             end

            ax = p.Results.ax;
            q = obj.(prop);

            plt = plot(ax, obj.xSurf(obj.xSurf>p.Results.xrange(1)&obj.xSurf<p.Results.xrange(2)), ...    % x vals withhin x range
                q(obj.xSurf>p.Results.xrange(1)&obj.xSurf<p.Results.xrange(2)), ...                       % corresponding y vals
                p.Results.fmt, ...                                                                      % Set format string
                'LineWidth', p.Results.LineWidth);                                                       % Set line width

            set(ax,'FontSize',p.Results.FontSize)

            xlim([x0 x1])
            ylim(p.Results.lims)                                                                        % Set y lims

%             if nargin > 5 && ~isempty(fmt)
%                 if isempty(xrange)
%                     plt = plot(ax,obj.xSurf,q,fmt);
%                 else
%                     plt = plot(ax,obj.xSurf(obj.xSurf>xrange(1)&obj.xSurf<xrange(2)),q(obj.xSurf>xrange(1)&obj.xSurf<xrange(2)),fmt);
%                 end
%             elseif nargin>4 && ~isempty(xrange)
%                 plt = plot(ax,obj.xSurf(obj.xSurf>xrange(1)&obj.xSurf<xrange(2)),q(obj.xSurf>xrange(1)&obj.xSurf<xrange(2)));
%             else
%                 plt = plot(ax,obj.xSurf,q);
%             end
%             if nargin > 3 && ~isempty(lims)
%                 ylim(lims);
%             end
%             set(plt,'LineWidth',1.5)
            disp('')
        end

        function value = thickness2ind(obj, del)
            y = obj.yO;

            if length(del) ~= size(y, 1)
                value = [];
                fprintf("Length of thickness vactor must equal lenghth of O grid\n");
            else
                for i=1:size(y,1)
                    [~, value(i)] = min(abs(y(i,:) - del(i)));
                end
            end

        end

        function value = edge_prop(obj, prop, varargin)

            p = inputParser;

            addParameter(p, 'thickness', []);

            parse(p, varargin{:});

            del = p.Results.thickness;

            if isempty(del)
                inds = obj.BLedgeInd;
            else
                if ischar(del)
                    del = obj.(del);
                end
    
                inds = obj.thickness2ind(del);
            end
            value = zeros(length(inds),1);
            q = obj.oGridProp(prop);
            for i=1:length(inds)
                value(i) = q(i, inds(i));
            end
        end

        function plt = plot_along_edge(obj, prop, varargin)

            x0 = min(obj.xSurf);
            x1 = max(obj.xSurf);

            defaultLineWidth = 1.5;
            defaultFmtString = '';
            defaultXRange = [x0 x1];

            p = inputParser;

            addParameter(p, 'ax', gca);
            addParameter(p, 'xrange', defaultXRange);
            addParameter(p, 'fmt', defaultFmtString);
            addParameter(p, 'LineWidth', defaultLineWidth);
            addParameter(p, 'ColorOrderIndex',[]);
            addParameter(p, 'thickness', 'delta99');

            parse(p, varargin{:});

            x = obj.xSurf;
            q = obj.edge_prop(prop, 'thickness', p.Results.thickness);
            i1 = obj.x2ind(p.Results.xrange(1));
            i2 = obj.x2ind(p.Results.xrange(2));

            plt = plot(p.Results.ax, x(i1:i2), q(i1:i2), p.Results.fmt, 'LineWidth',p.Results.LineWidth);

        end

        function plt = plot_BL_edge(obj, varargin)

            x0 = min(obj.xSurf);
            x1 = max(obj.xSurf);

            defaultLineWidth = 1.5;
            defaultFmtString = '';
            defaultXRange = [x0 x1];

            p = inputParser;

            addParameter(p, 'ax', gca);
            addParameter(p, 'xrange', defaultXRange);
            addParameter(p, 'fmt', defaultFmtString);
            addParameter(p, 'LineWidth', defaultLineWidth);
            addParameter(p, 'ColorOrderIndex',[]);
            addParameter(p, 'thickness', 'delta99');

            parse(p, varargin{:});

            x = obj.xSurf;
            if ischar(p.Results.thickness)
                edge = obj.(p.Results.thickness);
            else
                edge = p.Results.thickness;
            end

            i1 = obj.x2ind(p.Results.xrange(1));
            i2 = obj.x2ind(p.Results.xrange(2));

            plt = plot(p.Results.ax, x(i1:i2), edge(i1:i2), p.Results.fmt, 'LineWidth',p.Results.LineWidth);
            
        end

        function [f] = plot_BL_summary(obj, varargin)

            x0 = min(obj.xSurf);
            x1 = max(obj.xSurf);

            defaultLims = 'auto';
            defaultLineWidth = 1.5;
            defaultFmtString = '';
            defaultXRange = [x0-1 x1+1];
            if ~isempty(obj.label)
                defaultLabel = obj.label;
            else
                defaultLabel = [];
            end

            p = inputParser;

%             addRequired(p, 'prop');
            addParameter(p, 'fig', []);
            addParameter(p, 'lims', defaultLims);
            addParameter(p, 'xrange', defaultXRange);
            addParameter(p, 'loopxrange', []);
            addParameter(p, 'ploteq', false);
            addParameter(p, 'fmt', defaultFmtString);
            addParameter(p, 'LineWidth', defaultLineWidth);
            addParameter(p, 'label', defaultLabel);
            addParameter(p, 'ColorOrderIndex',[]);
            addParameter(p, 'FontSize', 12)
            addParameter(p, 'title', [])

            parse(p, varargin{:});

            if isempty(p.Results.fig)
                f = figure;
                t = tiledlayout(3,3,'TileSpacing','tight');
            else
                f = p.Results.fig;
                figure(f);
                if isempty(f.Children)
                    t = tiledlayout(3,3,'TileSpacing','tight');
                else
                    t = f.Children;
                end
            end

            if isempty(p.Results.loopxrange)
                loopxrange = p.Results.xrange;
            else
                loopxrange = p.Results.loopxrange;
            end
    
            try
                % Surface isentropic Mach number
                nexttile(1)
                a = gca;
                if ~isempty(p.Results.ColorOrderIndex)
                    a.ColorOrderIndex = p.Results.ColorOrderIndex;
                end
                hold on
                obj.blDevPlot('Msurf','xrange',p.Results.xrange,'fmt',p.Results.fmt,'LineWidth',p.Results.LineWidth, 'FontSize', p.Results.FontSize);
                ylabel('M_{surf}');
                a.XTickLabel = [];
                grid on
            catch
                disp('Could not plot Msurf')
            end
    
                try
                % Kinematic shape factor
                nexttile(2)
                a = gca;
                if ~isempty(p.Results.ColorOrderIndex)
                    a.ColorOrderIndex = p.Results.ColorOrderIndex;
                end
                hold on
                obj.blDevPlot('H_k','xrange',p.Results.xrange,'fmt',p.Results.fmt,'LineWidth',p.Results.LineWidth, 'FontSize', p.Results.FontSize);
                ylabel('H_k')
                a.XTickLabel = [];
                grid on
            catch
                disp('Could not plot Hk')
            end

            try
                % H-Pr loop
                nexttile(3,[2 1])
                hold on
                a = gca;
                if ~isempty(p.Results.ColorOrderIndex)
                    a = gca;
                    a.ColorOrderIndex = p.Results.ColorOrderIndex;
                end
                
    
                s = obj.plot_Hk_Pr_locus('xrange',loopxrange,'fmt',p.Results.fmt,'LineWidth',p.Results.LineWidth, 'FontSize', p.Results.FontSize);
                s.DisplayName = p.Results.label;
    
                if isempty(a.Legend)
                    l = legend(s,'Interpreter','latex');
                    l.Layout.Tile = 'east';
                end
    
                if p.Results.ploteq
                    lims = axis;
                    eq = obj.plot_Hk_Pr_locus('LineWidth',p.Results.LineWidth,'ploteq',true);
                    eq.DisplayName = 'Equilibrium Line';
                    axis(lims);
                    a.ColorOrderIndex = max(a.ColorOrderIndex - 1,1);
                    if length(a.Legend.Children) > 2 && isempty(find(strcmp(a.Legend.String, "Equilibrium Line"),1))
                        i = find(strcmp(a.Legend.String, "Equilibrium Line"));
                        a.Legend.Children = [a.Legend.Children(1:i-1) a.Legend.Children(i+1:end) a.Legend.Children(i)];
                    end
                end
            catch
                disp('Could not plot H-Pr loop')
            end

            try
                % Displacement thickness
                nexttile(4)
                a = gca;
                if ~isempty(p.Results.ColorOrderIndex)
                    a.ColorOrderIndex = p.Results.ColorOrderIndex;
                end
                    hold on
                obj.blDevPlot('delStar','xrange',p.Results.xrange,'fmt',p.Results.fmt,'LineWidth',p.Results.LineWidth, 'FontSize', p.Results.FontSize);
                ylabel('\delta^*/L');
                a.XTickLabel = [];
                grid on
            catch
                disp('Could not plot delStar')
            end

            try
                % Integrated production
                nexttile(5)
                a = gca;
                if ~isempty(p.Results.ColorOrderIndex)
                    a.ColorOrderIndex = p.Results.ColorOrderIndex;
                end
                hold on
                obj.blDevPlot('blPr','xrange',p.Results.xrange,'fmt',p.Results.fmt,'LineWidth',p.Results.LineWidth, 'FontSize', p.Results.FontSize);
                ylabel('Production')
                a.XTickLabel = [];
                grid on
            catch
                disp('Could not plot Pr')
            end

            try
                % Momentum thickness
                nexttile(7)
                if ~isempty(p.Results.ColorOrderIndex)
                    a = gca;
                    a.ColorOrderIndex = p.Results.ColorOrderIndex;
                end
                hold on
                obj.blDevPlot('theta','xrange',p.Results.xrange,'fmt',p.Results.fmt,'LineWidth',p.Results.LineWidth, 'FontSize', p.Results.FontSize);
                ylabel('\theta/L');
                xlabel('x/L');
                grid on
            catch
                disp('Could not plot theta')
            end

            try
                % Skin friction coeffient
                nexttile(8)
                a = gca;
                if ~isempty(p.Results.ColorOrderIndex)
                    a.ColorOrderIndex = p.Results.ColorOrderIndex;
                end
                hold on
                obj.blDevPlot('cf','xrange',p.Results.xrange,'fmt',p.Results.fmt,'LineWidth',p.Results.LineWidth, 'FontSize', p.Results.FontSize);
                ylabel('c_f')
                xlabel('x/L');
                grid on
            catch
                disp('Could not plot cf')
            end
            

            if ~isempty(p.Results.title)
                title(t, p.Results.title, 'FontSize', p.Results.FontSize+2);
            end

            nexttile(3)
            a = gca;
            a.Legend
            if isempty(a.Legend)
                l = legend(s,'Interpreter','latex');
                l.Layout.Tile = 'east';
            end

        end

        function plt = plot_y_profile(obj, x, prop, ax)
            if nargin < 4 || isempty(ax)
                ax = gca;
                disp('Creating axes')
            end

            [q, i] = BLprof(obj,x,prop);
            size(q);
            plot(ax, q, obj.yBL(i,:))
        end

        function s = plot_BL_profile(obj,x,prop, varargin)


            defaultAx = gca;
            defaultEdgeY = [];
            defaultLineWidth = 1.5;
            defaultFmtString = '';

            p = inputParser;

            addRequired(p, 'x', @isfloat);
            addRequired(p, 'prop', @ischar);
            addParameter(p, 'ax', defaultAx);
            addParameter(p, 'yEdge', defaultEdgeY);
            addParameter(p, 'normalise', false);
            addParameter(p, 'scale', 1);
            addParameter(p, 'fmt', defaultFmtString);
            addParameter(p, 'LineWidth', defaultLineWidth);
            addParameter(p, 'normaliseY', true);

            parse(p, x, prop, varargin{:});  

            [q, i] = BLprof(obj,x,prop);
            BLinds = obj.BLedgeInd;
            j = BLinds(i);
            if ~isempty(p.Results.yEdge)
                yEdge = p.Results.yEdge;
            elseif p.Results.normaliseY
                yEdge = obj.yBL(i,j);
            else
                yEdge = 1;
            end

            if p.Results.normalise
                scale = max(abs(q));
            else
                scale = p.Results.scale;
            end

            s = plot(p.Results.ax, q/scale, obj.yBL(i,:)/yEdge, ...                       % corresponding y vals
                p.Results.fmt, ...                                                                      % Set format string
                'LineWidth', p.Results.LineWidth);

            if ismember(string(prop),["dsdy","s"])
                hold on
                scatter(p.Results.ax, q(j)/scale, obj.yBL(i,j)/yEdge)
            end

        end

        function [plt] = plot_H_Pr_locus(obj, varargin) % ax, ploteq, xrange, fmt, lineColour)


            x0 = min(obj.xSurf);
            x1 = max(obj.xSurf);

            defaultAx = gca;
            defaultPlotEq = false;
            defaultXRange = [x0-1 x1+1];
            defaultFmtString = '';
            defaultLineWidth = 1.5;
            defaultLineColor = '';

            p = inputParser;

%             addRequired(p, 'prop');
            addParameter(p, 'ax', defaultAx);
            addParameter(p, 'ploteq', false)
            addParameter(p, 'xrange', defaultXRange);
            addParameter(p, 'fmt', defaultFmtString);
            addParameter(p, 'LineWidth', defaultLineWidth);
            addParameter(p, 'LineColor', '');

            parse(p, varargin{:})

            varplotargs = {};
            if p.Results.LineColor ~= ''
                varplotargs = [varplotargs {"Color" p.Results.LineColor}];
            end

            H = obj.H(obj.xSurf>p.Results.xrange(1)&obj.xSurf<p.Results.xrange(2));
            pr = obj.blPr(obj.xSurf>p.Results.xrange(1)&obj.xSurf<p.Results.xrange(2));

%             if nargin < 5
%                 locus_line = plot(ax,H,pr);
%             elseif nargin == 5 && ~isempty(fmt)
%                 locus_line = plot(ax,H,pr,fmt);
%             elseif nargin > 5 && ~isempty(lineColour) && ~isempty(fmt)
%                 fprintf('Colour specified\n')
%                 locus_line = plot(ax,H,pr,fmt,'Color',lineColour);
%             elseif nargin > 5 && ~isempty(lineColour) && isempty(fmt)
%                 locus_line = plot(ax,H,pr,'Color',lineColour);
%             end

            
            
            if p.Results.ploteq
                xtmp = linspace(1,3,51);% linspace(min(H),max(H),51);
                ytmp = 0.02456*((xtmp-1)./xtmp).^3;
                hold on
                plt = plot(xtmp,ytmp,'k:',"LineWidth", p.Results.LineWidth);
%                 legend([eq_line],'Equilibrium line','Location','northwest')
            else
                plt = plot(p.Results.ax, H, pr, "LineWidth", p.Results.LineWidth, varplotargs{:});
            end
            xlabel('H_{incomp}')
            ylabel('Pr')
            set(p.Results.ax,'FontSize',p.Results.FontSize)
            disp('')
            C = colororder;
            
        end

        function [plt] = plot_Hk_Pr_locus(obj, varargin) % ax, ploteq, xrange, fmt, lineColour)


            x0 = min(obj.xSurf);
            x1 = max(obj.xSurf);

            defaultAx = gca;
            defaultPlotEq = false;
            defaultXRange = [x0-1 x1+1];
            defaultFmtString = '';
            defaultLineWidth = 1.5;
            defaultLineColor = '';
            defaultFontSize = 12;

            p = inputParser;

%             addRequired(p, 'prop');
            addParameter(p, 'ax', defaultAx);
            addParameter(p, 'ploteq', false)
            addParameter(p, 'xrange', defaultXRange);
            addParameter(p, 'fmt', defaultFmtString);
            addParameter(p, 'LineWidth', defaultLineWidth);
            addParameter(p, 'LineColor', '');
            addParameter(p, 'FontSize', defaultFontSize);

            parse(p, varargin{:})

            fmt = p.Results.fmt;

%             varplotargs = {};
%             if p.Results.LineColor ~= ''
%                 varplotargs = [varplotargs {"Color" p.Results.LineColor}];
%             end

            H = obj.H_k(obj.xSurf>p.Results.xrange(1)&obj.xSurf<p.Results.xrange(2));
            pr = obj.blPr(obj.xSurf>p.Results.xrange(1)&obj.xSurf<p.Results.xrange(2));

%             if nargin < 5
%                 locus_line = plot(ax,H,pr);
%             elseif nargin == 5 && ~isempty(fmt)
%                 locus_line = plot(ax,H,pr,fmt);
%             elseif nargin > 5 && ~isempty(lineColour) && ~isempty(fmt)
%                 fprintf('Colour specified\n')
%                 locus_line = plot(ax,H,pr,fmt,'Color',lineColour);
%             elseif nargin > 5 && ~isempty(lineColour) && isempty(fmt)
%                 locus_line = plot(ax,H,pr,'Color',lineColour);
%             end

            
            
            if p.Results.ploteq
                xtmp = linspace(1,6,51);% linspace(min(H),max(H),51);
                ytmp = 0.02456*((xtmp-1)./xtmp).^3;
                hold on
                plt = plot(xtmp,ytmp,'k:',"LineWidth", p.Results.LineWidth);
%                 legend([eq_line],'Equilibrium line','Location','northwest')
            else
                plt = plot(p.Results.ax, H, pr, fmt, "LineWidth", p.Results.LineWidth);%varplotargs{:});
            end
            xlabel('H_{k}')
            ylabel('Pr')
            set(p.Results.ax,'FontSize',12)
            disp('')
            C = colororder;
            
        end

        function value = get.cPr(obj)
%             inds = obj.BLedgeInd;
%             Prnow = obj.oGridProp('Pr');
%             Unow = obj.U;
%             for i=1:size(obj.yBL,1)
%                 Prprof = Prnow(i,1:inds(i));
%                 Ue = Unow(i,inds(i));
%                 ys = obj.yBL(i,1:inds(i));
%                 value(i) = trapz(ys, Prprof)/Ue^3;
%             end

            inds = obj.BLedgeInd;
            Prnow = obj.oGridProp('Pr');
            ronow = obj.oGridProp('ro');
            Unow = obj.U;
            for i=1:size(obj.yBL,1)
                Ue = Unow(i,inds(i));
                roe = ronow(i, inds(i));
                ys = obj.yBL(i,1:inds(i));
                value(i,:) = Prnow(i,:)./(roe*Ue^3);
            end
            
        end

        function value = get.blPr(obj)
%             inds = obj.BLedgeInd;
%             Prnow = obj.oGridProp('Pr');
%             Unow = obj.U;
%             for i=1:size(obj.yBL,1)
%                 Prprof = Prnow(i,1:inds(i));
%                 Ue = Unow(i,inds(i));
%                 ys = obj.yBL(i,1:inds(i));
%                 value(i) = trapz(ys, Prprof)/Ue^3;
%             end

            inds = obj.BLedgeInd;
            Prnow = obj.oGridProp('Pr');
            ronow = obj.oGridProp('ro');
            Unow = obj.U;
            value(isnan(inds)) = NaN;
            for i=find(~isnan(inds)) %1:size(obj.yBL,1)
                if i==476
                    disp('')
                end
                Prprof = Prnow(i,1:inds(i));
                Ue = Unow(i,inds(i));
                roe = ronow(i, inds(i));
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, Prprof)/(roe*Ue^3);
            end
            value = obj.smooth_dist(value);
        end

        function value = get.blPr_dimensional(obj)
%             inds = obj.BLedgeInd;
%             Prnow = obj.oGridProp('Pr');
%             Unow = obj.U;
%             for i=1:size(obj.yBL,1)
%                 Prprof = Prnow(i,1:inds(i));
%                 Ue = Unow(i,inds(i));
%                 ys = obj.yBL(i,1:inds(i));
%                 value(i) = trapz(ys, Prprof)/Ue^3;
%             end

            inds = obj.BLedgeInd;
            Prnow = obj.oGridProp('Pr');
            ronow = obj.oGridProp('ro');
            Unow = obj.U;
            value(isnan(inds)) = NaN;
            for i=find(~isnan(inds)) %1:size(obj.yBL,1)
                if i==517
                    disp('')
                end
                Prprof = Prnow(i,1:inds(i));
                Ue = Unow(i,inds(i));
                roe = ronow(i, inds(i));
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, Prprof);%/(roe*Ue^3);
            end
            % value = obj.smooth_dist(value);
        end


        function value = get.cd_inner(obj)
            inds = obj.BLedgeInd;
            muSij2 = obj.oGridProp('nuSij2');
            Unow = obj.U;
            ronow = obj.oGridProp('ro');
            value(isnan(inds)) = NaN;
            for i=find(~isnan(inds))
                prof = muSij2(i,1:inds(i));
                Ue = Unow(i,inds(i));
                roe = ronow(i, inds(i));
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, prof)/(Ue^3);
            end
            value = obj.smooth_dist(value);
        end

        function value = get.cd_mean_strain(obj)
            inds = obj.BLedgeInd;
            diss_ave = obj.oGridProp('diss_ave');
            Unow = obj.U;
            ronow = obj.oGridProp('ro');
            value(isnan(inds)) = NaN;
            for i=find(~isnan(inds))
                prof = diss_ave(i,1:inds(i));
                Ue = Unow(i,inds(i));
                roe = ronow(i, inds(i));
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, prof)/(roe*Ue^3);
            end
            value = obj.smooth_dist(value);
        end

        function value = get.cd(obj)
            value = obj.cd_mean_strain + obj.blPr;
        end

        function value = get.blPr_eq(obj)
            Hk = obj.H_k;
            value = 0.02456*((Hk-1)./Hk).^3;
        end

        function value = get.Re_theta(obj)
            inds = obj.BLedgeInd;
            munow = obj.oGridProp('mu');
            
            ronow = obj.oGridProp('ro');
            size(ronow)
            size(inds)
            Unow = obj.U;
            thnow = obj.theta;
            value(isnan(inds)) = NaN;
            for i=find(~isnan(inds))
                value(i) = ronow(i,inds(i)).*Unow(i,inds(i)).*(thnow(i))./munow(i,inds(i));
            end

        end

        function value = get.Ue(obj)
            inds = obj.BLedgeInd;
            Unow = obj.U;                                              
            value(isnan(inds)) = NaN;
            for i=find(~isnan(inds))
                value(i) = Unow(i,inds(i));
            end
            value = obj.smooth_dist(value);
        end

        function value = get.Me(obj)
            inds = obj.BLedgeInd;
            Mnow = obj.oGridProp('M');                                             
            value(isnan(inds)) = NaN;
            for i=find(~isnan(inds))
                value(i) = Mnow(i,inds(i));
            end
            value = obj.smooth_dist(value);
        end

        function value = get.cf(obj)
            inds = obj.BLedgeInd;
            ronow = obj.oGridProp('ro');
            Unow = obj.U;
            roe(isnan(inds)) = NaN;
            Ue(isnan(inds)) = NaN;
            for i=find(~isnan(inds))
                roe(i) = ronow(i,inds(i));
                Ue(i) = Unow(i,inds(i));
            end
                value = obj.tau_w'./(0.5*roe.*Ue.*Ue);
        end


        
        function value = get.ctau(obj)
            value = obj.get_ctau;
        end

        % function value = get_ctau(obj)
        %     inds = obj.BLedgeInd;
        %     Uenow = obj.Ue;
        %     ronow = obj.oGridProp('ro');
        %     roe(isnan(inds)) = NaN;
        %     for i=find(~isnan(inds))
        %         roe(i) = ronow(i,inds(i));
        %     end
        % 
        %     Prnow = obj.oGridProp('Pr');
        %     dUdynow = obj.dUdy;
        %     tau = Prnow./dUdynow;
        %     % value = tau./(roe'.*Uenow.^2);
        %     value = tau./(ronow.*Uenow.^2);
        % end

        function value = get_ctau(obj)
            inds = obj.BLedgeInd;
            Uenow = obj.Ue;
            ronow = obj.oGridProp('ro');
            roe(isnan(inds)) = NaN;
            for i=find(~isnan(inds))
                roe(i) = ronow(i,inds(i));
            end
            
            Prnow = obj.oGridProp('Pr');
            dUdynow = obj.oGridProp('S_an_mag');
            tau = Prnow./dUdynow;
            % value = tau./(roe'.*Uenow.^2);
            value = tau./(ronow.*Uenow.^2);
        end

        function value = get.ct(obj)
            value = obj.blPr./(1-obj.Us);
        end

        function value = get.ctau_max(obj)
            ctau = obj.ctau;
            value = max(ctau,[],2);
        end
                
        function value = y_ctau_max(obj)
            ctau = obj.ctau;
            [~, inds] = max(ctau,[],2);
            value(isnan(inds)) = NaN;
            for i=find(~isnan(inds))
                value(i) = obj.yBL(i,inds(i));
            end
            value = obj.smooth_dist(value);
        end

        function value = j_ctau_max(obj)
            ctau = obj.ctau;
            [~, value] = max(ctau,[],2);
        end

                

        function value = get.pdyn(obj)
            inds = obj.BLedgeInd;
            ronow = obj.oGridProp('ro');
            Unow = obj.U;
            roe(isnan(inds)) = NaN;
            Ue(isnan(inds)) = NaN;
            for i=find(~isnan(inds))
                roe(i) = ronow(i,inds(i));
                Ue(i) = Unow(i,inds(i));
            end
            value = 0.5*roe.*Ue.*Ue;
        end

        function value = get.tau_w(obj)

            Unow = obj.U(:,2);
            Y0 = obj.yBL(:,2);
            munow = obj.oGridProp('mu');
            value = munow(:,2).*Unow./Y0;

        end

        function value = get.u_tau(obj)

            ronow = obj.oGridProp('ro');
            ro_w = ronow(:,1);

            value = sqrt(abs(obj.tau_w./ro_w));

        end

        function value = get.wall_scale(obj)
            nunow = obj.oGridProp('nu');
            nunow = nunow(:,2);
            value = nunow./obj.u_tau;
        end

        function [xplus,yplus,zplus] = wall_coords_offset(obj)
            dy = obj.yBL(:,2);
            size(dy)
            ds = obj.ssurf(2:end) - obj.ssurf(1:end-1);
            ds(end+1) = ds(end);
            ds = ds';
            size(ds)
            

            munow = obj.oGridProp('mu');
            munow = munow(:,2);
            ronow = obj.oGridProp('ro');
            ronow = ronow(:,2);

            xplus = ds.*sqrt(abs(obj.tau_w).*ronow)./munow;
            yplus = dy.*sqrt(abs(obj.tau_w).*ronow)./munow;
            if any(strcmp(properties(obj), 'span'))
                dz = ones(size(dy))*obj.span/(obj.nk-1);
                zplus = dz.*sqrt(abs(obj.tau_w).*ronow)./munow;
            else
                zplus = [];
            end
        end


        
        
        function blContour(obj, prop, varargin)

            defaultAx = gca;
            % defaultFmt = '';

            p = inputParser;

            addParameter(p, 'ax', defaultAx);
            addParameter(p, 'fmt', '');
            addParameter(p, 'scaleY', false);
            addParameter(p, 'lims', []);
            addParameter(p, 'label', prop);

            parse(p, varargin{:});

            ax = p.Results.ax;

            q = obj.oGridProp(prop);

            x = repmat(obj.xSurf, 1, size(obj.yBL,2));
            y = obj.yBL;
            if p.Results.scaleY
                y = y./obj.wall_scale;
            end

            pcolor(ax,x,y,q);
            shading interp
            if ~p.Results.scaleY
                axis equal
            end
            if ~isempty(p.Results.lims)
                caxis(p.Results.lims);
            end

            cb = colorbar(ax,"southoutside");
            cb.Label.String = p.Results.label;
            cb.Label.Interpreter = 'latex';
            cb.Label.FontSize = 12;

        end

        function [i, j, blk] = grid_inds_at_y_plus(obj,x,y_plus)
            [~,yplus_wall,~] = obj.wall_coords_offset;
            io = obj.x2ind(x);
            ys = obj.yBL(io,:);
            ynow = ys(2) * y_plus / yplus_wall(io);
            [~, jo] = min(abs(ys - ynow));
            i = obj.iO(io, jo);
            j = obj.jO(io, jo);
            blk = obj.blkO(io, jo);
        end

        function [i, j, blk] = grid_inds(obj,x,y)
            io = obj.x2ind(x);
            ys = obj.yBL(io,:);
            [~, jo] = min(abs(ys - y));
            i = obj.iO(io, jo);
            j = obj.jO(io, jo);
            blk = obj.blkO(io, jo);
        end

        function plotWallCoords(obj)

            [xplus,yplus,zplus] = obj.wall_coords_offset;
            yyaxis left
            plot(obj.xSurf,yplus);
            hold on
            yyaxis right
            plot(obj.xSurf,xplus);
            plot(obj.xSurf,zplus);
            xlabel('x/c')
            legend('y^+','x^+','z^+')

        end

        function [xplusm,yplusm,zplusm] = meanWallCoords(obj)
            [xplus,yplus,zplus] = obj.wall_coords_offset;
            xplusm = mean(xplus);
            yplusm = mean(yplus);
            zplusm = mean(zplus);
            fprintf('Mean x+: %4.2f\n', xplusm)
            fprintf('Mean y+: %4.2f\n', yplusm)
            fprintf('Mean z+: %4.2f\n', zplusm)
        end

        function plotYplus(obj)

            [~,yplus,~] = obj.wall_coords_offset;
            plot(obj.xSurf,yplus,'k-');
            xlabel('x/c')
            ylabel('y^+')
            grid on
            pbaspect([1 0.5 1])

        end

        function value = get.nu_e(obj)
            Tnow = obj.oGridProp('T');
            ronow = obj.oGridProp('ro');
            inds = obj.BLedgeInd;
            for i=1:size(Tnow,1)
                Te(i) = Tnow(i,inds(i));
                roe(i) = ronow(i,inds(i));
            end
            mu_e = obj.gas.mu_ref*(Te/obj.gas.mu_tref).^(3/2).*...
                (obj.gas.mu_tref+obj.gas.mu_cref)./(Te+obj.gas.mu_cref);
            value = mu_e./roe;
            value = obj.smooth_dist(value);
        end

        function value = mass_average_inlet(obj, prop, i)

            blks = obj.blk.inlet_blocks{1};

            if nargin < 3
                i = 20;
            end
            
            flux = 0;
            mass = 0;

            for ib = blks

                xnow = obj.blk.x{ib}(i,:);
                ynow = obj.blk.y{ib}(i,:);

                rounow = obj.ro{ib}(i,:).*obj.u{ib}(i,:);
                rovnow = obj.ro{ib}(i,:).*obj.v{ib}(i,:);

                xfluxnow = obj.ro{ib}(i,:).*obj.u{ib}(i,:).*obj.(prop){ib}(i,:);
                yfluxnow = obj.ro{ib}(i,:).*obj.v{ib}(i,:).*obj.(prop){ib}(i,:);

                mass = mass + trapz(ynow,rounow) - trapz(xnow, rovnow);
                flux = flux + trapz(ynow, xfluxnow) - trapz(xnow, yfluxnow);

            end

            value = flux/mass;
            
        end

        function value = mass_average_outlet(obj, prop, i)

            blks = obj.blk.outlet_blocks{1};

            if nargin < 3
                i = obj.blk.blockdims(blks(1), 1) - 20;
            else 
                i = obj.blk.blockdims(blks(1), 1) - i;
            end
            
            flux = 0;
            mass = 0;

            for ib = blks

                xnow = obj.blk.x{ib}(i,:);
                ynow = obj.blk.y{ib}(i,:);

                rounow = obj.ro{ib}(i,:).*obj.u{ib}(i,:);
                rovnow = obj.ro{ib}(i,:).*obj.v{ib}(i,:);

                xfluxnow = obj.ro{ib}(i,:).*obj.u{ib}(i,:).*obj.(prop){ib}(i,:);
                yfluxnow = obj.ro{ib}(i,:).*obj.v{ib}(i,:).*obj.(prop){ib}(i,:);

                mass = mass + trapz(ynow,rounow) - trapz(xnow, rovnow);
                flux = flux + trapz(ynow, xfluxnow) - trapz(xnow, yfluxnow);

            end

            value = flux/mass;
            
        end

        function value = area_average_inlet(obj, prop, i)

            blks = obj.blk.inlet_blocks{1};

            if nargin < 3
                i = 20;
            end
            
            flux = 0;
            area = 0;

            for ib = blks

                xnow = obj.blk.x{ib}(i,:);
                ynow = obj.blk.y{ib}(i,:);

                dx = diff(xnow);
                dy = diff(ynow);
                ds = sqrt(dx.^2+dy.^2);
                s = [0 cumsum(ds)];
                
                area = area + trapz(s, ones(size(s)));
                flux = flux + trapz(s, obj.(prop){ib}(i,:));

            end

            value = flux/area;
            
        end
        
        function value = area_average_outlet(obj, prop, i)

            blks = obj.blk.outlet_blocks{1};

            if nargin < 3
                i = obj.blk.blockdims(blks(1), 1) - 20;
            end
            
            flux = 0;
            area = 0;

            for ib = blks

                xnow = obj.blk.x{ib}(i,:);
                ynow = obj.blk.y{ib}(i,:);

                dx = diff(xnow);
                dy = diff(ynow);
                ds = sqrt(dx.^2+dy.^2);
                s = [0 cumsum(ds)];

                area = area + trapz(s, 1);
                flux = flux + trapz(s, obj.(prop){ib}(i,:));

            end

            value = flux/area;
            
        end

        function value = get.alpha1(obj)
            inlet_blks = obj.blk.inlet_blocks{1};
            xmom = 0;
            ymom = 0;
            alp = [];
            ynow = [];
            ronow = [];
            den = 0;
            num = 0;
            mass = 0;
            roave = 0;
            for ib = inlet_blks
                i=20;%size(obj.blk.x{ib},1);
%                 ynow = [ynow obj.blk.y{ib}(i,:)];
%                 alp = [alp obj.v{ib}(i,:)./obj.u{ib}(i,:)];
%                 ronow = [ronow obj.ro{ib}(i,:)];
                ynow = obj.blk.y{ib}(i,:);
%                 anow = atand(obj.u{ib}(i,:)./obj.v{ib}(i,:));
                rounow = obj.ro{ib}(i,:).*obj.u{ib}(i,:);
                rouvnow = obj.ro{ib}(i,:).*obj.u{ib}(i,:).*obj.v{ib}(i,:);
%                 mnow = obj.ro{ib}(i,:).*obj.u{ib}(i,:);
%                 num = num+trapz(ynow,anow.*mnow);
%                 den = den+trapz(ynow,mnow);
                mass = mass + trapz(ynow,rounow);
                roave = roave + trapz(ynow, obj.ro{ib}(i,:));
                ymom = ymom + trapz(ynow,rouvnow);
            end
%             [ynow, is] = sort(ynow);
%             ronow = ronow(is);
%             alp = alp(is);
            v = ymom/mass;
            u = mass/roave;

            value = atan2d(v,u);
%             value = num/den;
            

        end

        function value = get.alpha2(obj)
            outlet_blks = obj.blk.outlet_blocks{1};
            xmom = 0;
            ymom = 0;
            alp = [];
            ynow = [];
            ronow = [];
            den = 0;
            num = 0;
            mass = 0;
            roave = 0;
            for ib = outlet_blks
                i=size(obj.blk.x{ib},1);
%                 ynow = [ynow obj.blk.y{ib}(i,:)];
%                 alp = [alp obj.v{ib}(i,:)./obj.u{ib}(i,:)];
%                 ronow = [ronow obj.ro{ib}(i,:)];
                ynow = obj.blk.y{ib}(i,:);
%                 anow = atand(obj.u{ib}(i,:)./obj.v{ib}(i,:));
                rounow = obj.ro{ib}(i,:).*obj.u{ib}(i,:);
                rouvnow = obj.ro{ib}(i,:).*obj.u{ib}(i,:).*obj.v{ib}(i,:);
%                 mnow = obj.ro{ib}(i,:).*obj.u{ib}(i,:);
%                 num = num+trapz(ynow,anow.*mnow);
%                 den = den+trapz(ynow,mnow);
                mass = mass + trapz(ynow,rounow);
                roave = roave + trapz(ynow, obj.ro{ib}(i,:));
                ymom = ymom + trapz(ynow,rouvnow);
            end
%             [ynow, is] = sort(ynow);
%             ronow = ronow(is);
%             alp = alp(is);
            v = ymom/mass;
            u = mass/roave;

            value = atan2d(v,u);
%             value = num/den;
            

        end

        function value = get.Pr_nondim(obj)

            inds = obj.BLedgeInd;
            Prnow = obj.oGridProp('Pr');
            ronow = obj.oGridProp('ro');
            Unow = obj.oGridProp('U');
            delnow = obj.delta99;

            for i=1:length(inds)
                value(i,:) = Prnow(i,:)*delnow(i)/(ronow(i, inds(i)) * Unow(i, inds(i))^3);
            end

        end

        function p = plot_wake(obj, prop, x, fmt)

	    if nargin < 4
		    fmt = '';
	    end
            i = find_nearest(obj.blk.x{obj.blk.outlet_blocks{1}(1)}(:,1), x);
    
            y = [];
            prof = [];

            switch prop
                case 'zeta'

                    h0in = obj.mass_average_inlet('h0', 20);
                    hin = obj.mass_average_inlet('h', 20);
                    sin = obj.mass_average_inlet('s', 50);
                    T2 = obj.mass_average_outlet('T', i);
                    scale = T2/(h0in-hin);

                    tmp = obj.s;
                    offset = sin;
                case 'alpha'
                    tmp = obj.flowAng;
                    scale=1;
                    offset=0;
	    	case 'alpha_norm'
		    tmp = obj.flowAng
		    scale = 1;
		    offset = obj.alpha2;
	    
            end
        
        
            for ib = obj.blk.outlet_blocks{1}
                y = [y obj.blk.y{ib}(i,1:end-1)];
                prof = [prof tmp{ib}(i,1:end-1)];
            end

            prof = scale*(prof-offset);
        
            [yprof, inds] = sort(y);
            prof = prof(inds);

            yprof = (yprof-yprof(1))/(yprof(end)-yprof(1));

            p = plot(yprof, prof, fmt)

        end



        function value = get.zeta(obj)
            Tnow = obj.mass_average_outlet('T', 20);
            dsnow = obj.mass_average_outlet('s', 20) - obj.mass_average_inlet('s', 20);
            h0now = obj.mass_average_inlet('h0', 20);
            hnow = obj.mass_average_inlet('h', 20);

            value = (Tnow*dsnow)/(h0now-hnow);
                    
        end

        function value = pout(obj)
            outlet_blks = obj.blk.outlet_blocks{1};
            xmom = 0;
            ymom = 0;
            alp = [];
            ynow = [];
            ronow = [];
            den = 0;
            num = 0;
            mass = 0;
            roave = 0;
            roup = 0;
            for ib = outlet_blks
                i=size(obj.blk.x{ib},1);
%                 ynow = [ynow obj.blk.y{ib}(i,:)];
%                 alp = [alp obj.v{ib}(i,:)./obj.u{ib}(i,:)];
%                 ronow = [ronow obj.ro{ib}(i,:)];
                ynow = obj.blk.y{ib}(i,:);
%                 anow = atand(obj.u{ib}(i,:)./obj.v{ib}(i,:));
                rounow = obj.ro{ib}(i,:).*obj.u{ib}(i,:);
                roupnow = obj.ro{ib}(i,:).*obj.u{ib}(i,:).*obj.p{ib}(i,:);
                % rouvnow = obj.ro{ib}(i,:).*obj.u{ib}(i,:).*obj.v{ib}(i,:);
%                 mnow = obj.ro{ib}(i,:).*obj.u{ib}(i,:);
%                 num = num+trapz(ynow,anow.*mnow);
%                 den = den+trapz(ynow,mnow);
                mass = mass + trapz(ynow,rounow);
                roup = roup + trapz(ynow, roupnow);
                % ymom = ymom + trapz(ynow,rouvnow);
            end
%             [ynow, is] = sort(ynow);
%             ronow = ronow(is);
%             alp = alp(is);
            % v = ymom/mass;
            % u = mass/roave;

            value = roup/mass;
%             value = num/den;
            

        end

        function value = mdot(obj, bnd)
            if nargin < 2
                bnd = 'outlet';
            end

            switch bnd
                case 'inlet'
                    blks = obj.blk.inlet_blocks{1};
                    i=20;
                case 'outlet'
                    blks = obj.blk.outlet_blocks{1};
                    i = obj.blk.blockdims(blks(1), 1) - 20;
            end

            value = 0;

            for ib=blks

                xnow = obj.blk.x{ib}(i,:);
                ynow = obj.blk.y{ib}(i,:);

                ru = obj.ro{ib}(i,:).*obj.u{ib}(i,:);
                rv = obj.ro{ib}(i,:).*obj.v{ib}(i,:);
                
                % dx = xnow(2:end)-xnow(1:end-1);
                % dy = ynow(2:end)-ynow(1:end-1);
                % ru = 0.5 * (runow(2:end) + runow(1:end-1));
                % rv = 0.5 * (rvnow(2:end) + rvnow(1:end-1));
                
                value = value + trapz(ynow, ru) - trapz(xnow, rv);
            end

        end

        function value = get.p0out(obj)
            
            outlet_blks = obj.blk.outlet_blocks{1};

            num = 0;
            mass = 0;
            for ib = outlet_blks
                i=ceil(0.95*size(obj.blk.x{ib},1));
                ynow = obj.blk.y{ib}(i,:);
                p0now = obj.p0{ib}(i,:);
                rounow = obj.ro{ib}(i,:).*obj.u{ib}(i,:);
                mass = mass + trapz(ynow,rounow);
                num = num + trapz(ynow, rounow.*p0now);
            end

            value = num/den;

        end

        function value = get.yplus(obj)
                        dy = obj.yBL(:,2);
            size(dy)
            ds = obj.ssurf(2:end) - obj.ssurf(1:end-1);
            ds(end+1) = ds(end);
            ds = ds';
            size(ds)
            

            munow = obj.oGridProp('mu');
            munow = munow(:,2);
            ronow = obj.oGridProp('ro');
            ronow = ronow(:,2);

            value = obj.yBL.*sqrt(abs(obj.tau_w).*ronow)./munow;

        end

        function value = get.mut(obj)
            value = obj.get_mut;
        end

        function obj = set.mut(obj, value)
            obj.set_mut(value);
        end

        function set_mut(obj,value)
        end

        function value = get_mut(obj)
            disp('Overload get_mut in relevant subclass')
            value = [];
        end

        function value = get_mut_ratio(obj)
            mutnow = obj.mut;
            munow = obj.mu;
            for ib = 1:obj.NB
                value{ib} = mutnow{ib}./munow{ib};
            end
        end

        function value = get.mut_ratio(obj)
            value = obj.get_mut_ratio;
        end

        function value = get.Yp(obj)
            p0now = obj.p0;
            
            for ib=1:obj.NB
                value{ib} = (obj.p0in - p0now{ib})/(obj.p0in - obj.pin);
            end

        end

        function value = bl_history(obj, xrange)

            if nargin < 2 || isempty(xrange)
                xrange = [obj.xSurf(1) obj.xSurf(end)];
            end

            is = obj.x2ind(xrange(1)):obj.x2ind(xrange(2));

            props = ["H", "H_k", "H_ke", "H_rho", "delStar", "theta",...
                "Re_theta", "cf", "cd", "ctau_max", "Ue", "Us", ...
                "ct", "blPr", "xSurf", "Msurf", "ctau", "wall_scale", ...
                "yBL", "delta"];

            for i=1:length(props)
                q = obj.(props(i));
                if size(q,1) == 1
                    value.(props(i)) = q(is);
                else
                    value.(props(i)) = q(is,:);
                end
            end
            value.x = value.xSurf;
            value.Pr = value.blPr;
            value.Pk = obj.oGridProp('Pr');
            value.roe = obj.edge_prop('ro');
            value.bl_inds = obj.BLedgeInd;

        end

        function terms = shear_lag_terms(obj)

            

        end

        function value = get.pr_aligned_s(obj)
            tau_aligned = obj.align_tensor_with_surface('tau_Re');
            St_aligned = obj.align_tensor_with_surface('St');
            value = squeeze(tau_aligned(:,:,1,1)).*squeeze(St_aligned(:,:,1,1));
        end

        function value = get.pr_aligned_n(obj)
            tau_aligned = obj.align_tensor_with_surface('tau_Re');
            St_aligned = obj.align_tensor_with_surface('St');
            value = squeeze(tau_aligned(:,:,2,2)).*squeeze(St_aligned(:,:,2,2));
        end

        function value = get.pr_aligned_t(obj)
            tau_aligned = obj.align_tensor_with_surface('tau_Re');
            St_aligned = obj.align_tensor_with_surface('St');
            value = squeeze(tau_aligned(:,:,1,2)).*squeeze(St_aligned(:,:,1,2)) ...
            +  squeeze(tau_aligned(:,:,2,1)).*squeeze(St_aligned(:,:,2,1));
        end

        function sol = run_mrchbl(obj, xstart, path, xtrip, xmatch)

            Kcorr = 5.6;
            Kp = 1.0;
            Kd  = 1.0;

            if nargin < 5
                xmatch = [];
            end
            
            ni = 401;
            Min = obj.x2prop(xstart, 'Msurf');
            istart = obj.x2ind(xstart);

            rt_dns = obj.Re_theta(istart:end);
            ds_dns = obj.delStar(istart:end);
            th_dns = obj.theta(istart:end);

            rt = rt_dns(1);
            ds = ds_dns(1);
            th = th_dns(1);

            % ds = obj.x2prop(xstart, 'delStar');
            % th = obj.x2prop(xstart, 'theta');
            a0 = sqrt(obj.gas.gam*obj.gas.rgas*obj.bcs.Toin);
            % Ue = obj.Ue(istart:end)/a0;
            x = obj.xSurf(istart:end);
            xc = linspace(x(1), x(end), ni);
            xnow = linspace(x(1), x(end), floor(ni/4)+1);
            Msurf = obj.Msurf(istart:end);
            Mc = interp1(x, Msurf, xnow);
            M = interp1(xnow, Mc, xc, 'pchip');
            % Uc = interp1(x, Ue, xnow);
            % Ue = interp1(xnow, Uc, xc, 'pchip');

            Ue = V_MT0(M, obj.bcs.Toin)/a0;
            
            if nargin < 4 || isempty(xtrip)
                xtrip = 0.5*(xc(1)+xc(2));
            end
            if xc(1) == 0
                xc(1) = 0.1*xc(2);
            end

            H = ds/th;
            Hk = MISES_correlations.fHk(Min, H);
            Hks = MISES_correlations.fHks(Hk, rt);
            Hs = MISES_correlations.fHs(Min, Hks);
            Pr = obj.x2prop(xstart, 'blPr');

            Us = MISES_correlations.fUs(Hs, Hk, H);
            ct = Pr/(1-Us);

            input_file = fullfile(path, 'inp.dat');
            output_file = fullfile(path, 'out.dat');
            write_mrchbl_input(input_file, xc, Ue, Min, rt, xtrip, th, ds, ct, Kcorr, Kp, Kd);
            mrchbl_path = '~/MISES/MISES/bin/mblrun';
            system([mrchbl_path ' ' input_file ' ' output_file]);
            sol = read_mrchbl_output(output_file);

            if ~isempty(xmatch)

                dst = interp1(x, ds_dns, xmatch);
                tht = interp1(x, th_dns, xmatch);
                rtt = interp1(x, rt_dns, xmatch);

                for i=1:5
                    dsm = interp1(sol.x, sol.delStar, xmatch);
                    thm = interp1(sol.x, sol.theta, xmatch);
                    rtm = interp1(sol.x, sol.Re_theta, xmatch);
    
                    ds = ds - dsm + dst;
                    th = th - thm + tht;
                    rt = rt * (rtt/rtm);
    
                    write_mrchbl_input(input_file, xc, Ue, Min, rt, xtrip, th, ds, ct, Kcorr, Kp, Kd);
                    system([mrchbl_path ' ' input_file ' ' output_file]);
                    sol = read_mrchbl_output(output_file);
                end
            end

            
        end

        function sol = MISES_BL(obj, xstart, varargin)

            data = [];
            istart = obj.x2ind(xstart);
            xsurf = obj.xSurf;
            is = istart:length(xsurf);
            xsurf = xsurf(is);
            nue = obj.nu_e(is);

            Kdefault = struct('Kcorr', 4.2, 'Kd', 0, 'Kp', 0);

            p = inputParser;
            addParameter(p, 'K', Kdefault);
            addParameter(p, 'Me', []);
            addParameter(p, 'Ue', []);
            addParameter(p, 'DNSinputs',[]);

            parse(p, varargin{:});

            Me = p.Results.Me;
            Ue = p.Results.Ue;
            
            data.Kcorr = p.Results.K.Kcorr;
            data.Kd = p.Results.K.Kd;
            data.Kp = p.Results.K.Kp;

            v = string(p.Results.DNSinputs);

            if isempty(Me)
                Me = obj.Msurf(is);
            else
                Me = Me(is);
            end

            if isempty(Ue)
                Ue = obj.Ue(is);
            else
                Ue = Ue(is);
            end

            if ismember("theta", v) || ismember("thick", v)
                data.theta = obj.theta(is);
                theta0 = data.theta(1);
            else
                theta0 = obj.x2prop(xstart, 'theta');
            end

            Reth0 = Ue(1)*theta0/nue(1);

            if ismember("delStar", v) || ismember("thick", v)
                data.delStar = obj.delStar(is);
                delStar0 = data.delStat(1);
            else
                delStar0 = obj.x2prop(xstart, 'delStar');
            end

            if ismember("H", v) || ismember("SFs", v)
                data.H = obj.H(is);
            end
            H0 = delStar0/theta0;

            if ismember("Hk", v) || ismember("SFs", v)
                data.Hk = obj.H_k(is);
                Hk0 = data.Hk(1);
            else
                Hk0 = MISES_correlations.fHk(Me(1), H0);
            end

            Hks0 = MISES_correlations.fHks(Hk0, Reth0);

            if ismember("H_ke", v) || ismember("SFs", v)
                data.Hs = obj.H_ke(is);
                Hs0 = data.Hs(1);
            else
                Hs0 = MISES_correlations.fHs(Me(1), Hks0);
            end
            
            if ismember("cf", v)
                data.cf = obj.cf(is);
            end


            fun = @(x, y) MISES_ydash(x, y, xsurf, Me, Ue, nue, data);

            Pr0 = obj.x2prop(xstart, 'blPr');
            Us0 = MISES_correlations.fUs(Hs0, Hk0, H0);

            y0 = [theta0; ...
                Hs0; ...
                Pr0/(1-Us0)];

            step = xsurf(2)-xsurf(1);
            opts = odeset('InitialStep',step, 'MaxStep', 2*step, 'Stats','on');
            sol = ode78(fun, xsurf, y0, opts);

            sol.Me = interp1(xsurf, Me, sol.x);
            sol.Ue = interp1(xsurf, Ue, sol.x);
            sol.nue = interp1(xsurf, nue, sol.x);

            if ismember("theta", v) || ismember("thick", v)
                sol.theta = interp1(xsurf, data.theta, sol.x);
            else
                sol.theta = sol.y(1,:);
            end

            sol.Reth = sol.theta.*sol.Ue./sol.nue;

            if ismember("Hs", v) || ismember("SFs", v)
                sol.Hs = interp1(xsurf, data.Hs, sol.x);
            else
                sol.Hs = sol.y(2,:);
            end

            sol.Hks = MISES_correlations.fHks_Hs(sol.Hs, sol.Me);

            if ismember("Hk", v) || ismember("SFs", v)
                sol.H_k = interp1(xsurf, data.Hk, sol.x);
            else
                sol.H_k = MISES_correlations.fHk_Hks(sol.Hks, sol.Reth);
            end

            if ismember("H", v)
                sol.H = interp1(xsurf, data.H, sol.x);
            else
                sol.H = MISES_correlations.fH(sol.Me, sol.H_k);
                % sol.H = sol.delStar/sol.theta;
            end

            if ismember("delStar", v) || ismember("thick", v)
                sol.delStar = interp1(xsurf, data.delStar, sol.x);
            else
                sol.delStar = sol.H.*sol.theta;
            end

            sol.H_ks = MISES_correlations.fHks(sol.H_k, sol.Reth);
            

            if ismember("H_ke", v) || ismember("SFs", v)
                sol.H_ke = interp1(xsurf, data.Hs, sol.x);
            else
                sol.H_ke = sol.Hs; %MISES_correlations.fHs(sol.Me, sol.H_ks);
            end

            sol.Hss = MISES_correlations.fHss(sol.Me, sol.H_k);
            sol.ctau = sol.y(3,:);

            if ismember("cf", v)
                sol.cf = interp1(xsurf, data.cf, sol.x);
            else
                sol.cf = MISES_correlations.fCf(sol.H_k, sol.Reth, sol.Me);
            end

            sol.Us = MISES_correlations.fUs(sol.H_ke, sol.H_k, sol.H);
            sol.cd = MISES_correlations.fCd(sol.cf, sol.Us, sol.ctau);
            sol.CtauEq = MISES_correlations.fCtEQ(sol.H_ke, sol.Us, sol.H, sol.H_k);
            sol.del = MISES_correlations.fDel(sol.theta, sol.H_k, sol.delStar);
            sol.Pr = sol.ctau.*(1-sol.Us);
            sol.PrEq = sol.CtauEq.*(1-sol.Us);

        end        
    end
end
