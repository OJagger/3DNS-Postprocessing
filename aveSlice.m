classdef aveSlice < kCut
    % MEANSLICE Contains the 2D (spanwise averaged) mean flow

        
    properties
        p0in;
        Uinf;
        roinf;
        muinf;
        Tinf;
        T0in;
        dsdyThresh = -50;
        use_unsflo = false;
        blEdgeMode;
        iSmoothBL =  false;
        nSmooth;
    end

    properties (Dependent = true)
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
        Res;            % Surface distance Reynolds No
        blPr;           % Componant of cd due to production of tke
        blPr_eq;        % Coles' equilibrium prodution
        tau_w;          % Wall shear stress
        u_tau;          % Friction velocity
        dUdy;
        dTdy;
        cf;
        ctau;
        ctau_max;
        Re_theta;
        iPS;
        iMaxPr;
        iEq;
        Re_theta_ps;
        pdyn            % Dynamic pressure, 0.5*rho*V^2
        nu_e            % Boundary layer edge viscosity
        alpha2;
        p0out;
        yplus;
    end

    methods
        function obj = aveSlice(blk, gas, bcs)
            obj@kCut(blk, gas, bcs);
            disp('Constructing aveSlice')
            obj.nSmooth = 1;
        end

        function getBCs(obj, inlet_blocks, is)
            if nargin < 3
                is = 40:100;
            end
            if nargin < 2
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
            end
            obj.p0in = mean(p0,'all');
            obj.Tinf = mean(Tinf, 'all');
            obj.Uinf = mean(Uinf,'all');
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
       

        function value = BLedgeInd(obj, mode)
            if nargin < 2 || isempty(mode)
                mode = "sThresh";
            end
            u = obj.oGridProp('U');
            fprintf('BL edge detection mode: %s\n',mode)
            switch mode
                case "sGradThresh"
                    temp = obj.dsdy;
                    for i=1:size(temp,1)
                        j=size(temp,2);
                        while temp(i,j) > obj.dsdyThresh && j>1
                            j = j-1;
                        end
                        value(i) = min(j+1,size(temp,2));
                        [~, value(i)] = max(u(i,1:value(i)));
                    end
                case "sThresh"
                    temp = obj.oGridProp('s');
                    for i=1:size(temp,1)
                        j=size(temp,2);
                        se = temp(i,end);
                        sw = temp(i,1);
                        se;
                        while temp(i,j) < se + 0.02*(sw-se) && j>1
                            j = j-1;
                        end
                        value(i) = min(j+1,size(temp,2));
                        %[~, value(i)] = max(u(i,1:value(i)));
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
            end
            if obj.iSmoothBL
                value = obj.smooth_bl_edge(obj.xSurf,value);
            end
        end

        function [xedge, yedge] = getBLedgeCoords(obj,mode)
            if nargin < 2
                mode = "sGradThresh";
            end
            inds = obj.BLedgeInd(mode);
            for i=1:length(inds)
                xedge(i) = obj.xO(i,inds(i));
                yedge(i) = obj.yO(i,inds(i));
            end
        end

        function value = get.iPS(obj)
            x = obj.xSurf;
            M = smooth(obj.Msurf, obj.nSmooth);
            pr = smooth(obj.blPr, obj.nSmooth);

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
            pr = smooth(obj.blPr, obj.nSmooth);
            while pr(i) > pr(i-1)
                i = i+1;
            end
            value = i;
        end

        function value = delta99_unsflo(obj)
            %inds = obj.BLedgeInd;
            
        end

        function value = delta99(obj, mode)
            if nargin < 2
                mode = "sThresh";
            end
            switch mode
                case "sGradThresh"
                    inds = obj.BLedgeInd("sGradThresh");
                    for i=1:size(obj.yBL,1)
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
                    for i=1:size(obj.yBL,1)
                        value(i) = obj.yBL(i,inds(i));
                    end
            end
            value = smooth(value,obj.nSmooth);
            
        end


        function value = get.U(obj)
            disp('Calculating U')
            unow = obj.oGridProp('u');
            vnow = obj.oGridProp('v');
            nnow = obj.n;
            value = zeros(size(obj.yBL));
            for i=1:size(obj.yBL,1)
                tang = [-nnow(2,i); nnow(1,i)];
                for j=1:size(obj.yBL,2)
                    velnow = [unow(i,j); vnow(i,j)];
                    % value(i,j) = -dot(tang, velnow - nnow(:,i)*dot(nnow(:,i),velnow));
                    value(i,j) = norm(velnow - dot(velnow,nnow(:,i)));
                end
            end
        end

        function value = get.delStar(obj)
            inds = obj.BLedgeInd(obj.blEdgeMode);
            ronow = obj.oGridProp('ro');
            Unow = obj.U;
            value = zeros(1,length(inds));
            for i=1:size(obj.yBL,1)
                integrand = 1 - ronow(i,1:inds(i)).*Unow(i,1:inds(i))./(ronow(i,inds(i))*Unow(i,inds(i)));
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, integrand);
            end
            value = smooth(value,obj.nSmooth);
        end

        function value = get.delStar_k(obj)
            inds = obj.BLedgeInd;
            Unow = obj.U;
            value = zeros(1,length(inds));
            for i=1:size(obj.yBL,1)
                integrand = 1 - Unow(i,1:inds(i))./Unow(i,inds(i));
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, integrand);
            end
            value = smooth(value,obj.nSmooth);
        end

        function value = get.delRho(obj)
            inds = obj.BLedgeInd;
            ronow = obj.oGridProp('ro');
            Unow = obj.U;
            for i=1:size(obj.yBL,1)
                roprof = ronow(i,1:inds(i));
                Uprof = Unow(i,1:inds(i));
                ro0 = ronow(i,inds(i));
                U0 = Unow(i,inds(i));
                integrand = (1 - roprof/ro0).*Uprof/U0;
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, integrand);
            end
            value = smooth(value,obj.nSmooth);
        end

        function value = get.theta(obj)
            inds = obj.BLedgeInd;
            ronow = obj.oGridProp('ro');
            Unow = obj.U;
            for i=1:size(obj.yBL,1)
                roprof = ronow(i,1:inds(i));
                Uprof = Unow(i,1:inds(i));
                ro0 = ronow(i,inds(i));
                U0 = Unow(i,inds(i));
                integrand = (roprof.*Uprof/(ro0*U0)).*(1-Uprof/U0);
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, integrand);
            end
            value = smooth(value,obj.nSmooth);
        end

        function value = get.theta_k(obj)
            inds = obj.BLedgeInd;
            Unow = obj.U;
            for i=1:size(obj.yBL,1)
                Uprof = Unow(i,1:inds(i));
                U0 = Unow(i,inds(i));
                integrand = (Uprof/U0).*(1-Uprof/U0);
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, integrand);
            end
            value = smooth(value,obj.nSmooth);
        end

        function value = get.thetaStar(obj)
            inds = obj.BLedgeInd;
            ronow = obj.oGridProp('ro');
            Unow = obj.U;
            for i=1:size(obj.yBL,1)
                roprof = ronow(i,1:inds(i));
                Uprof = Unow(i,1:inds(i));
                ro0 = ronow(i,inds(i));
                U0 = Unow(i,inds(i));
                integrand = (roprof.*Uprof/(ro0*U0)).*(1-Uprof.^2/U0^2);
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, integrand);
            end
            value = smooth(value,obj.nSmooth);
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
            for i=1:size(obj.yBL,1)
                roprof = ronow(i,1:inds(i));
                Uprof = Unow(i,1:inds(i));
                ro0 = ronow(i,inds(i));
                U0 = Unow(i,inds(i));
                integrand = (roprof.*Uprop)/(ro0*U0);
                ys = obj.yBL(i,1:inds(i));
                num(i) = trapz(ys, integrand);
            end
            value = num(i)./obj.theta;
        end

        function value = get.delta(obj)
            value = (obj.H + obj.H1).*obj.theta;
        end

        function value = get.Res(obj)
            value = obj.ssurf*obj.Uinf*obj.roinf/obj.muinf;
        end

        function plot_BL_profile(obj,x,prop,ax)
            if nargin < 4 || isempty(ax)
                ax = gca;
                disp('Creating axes')
            end

            [q, i] = BLprof(obj,x,prop);
            size(q);
            q = q/max(q);
            BLinds = obj.BLedgeInd;
            if string(prop) == "dsdy"
                plot(ax, q, obj.yBL(i,:))
                hold on
                scatter(ax, q(BLinds(i)), obj.yBL(i,BLinds(i)))
            else
                plot(ax, q, obj.yBL(i,:)/obj.yBL(i,end))
                %hold on
                %scatter(ax, q(BLinds(i)), obj.yBL(i,BLinds(i))/obj.yBL(i,end))
            end
        end

        function plt = blDevPlot(obj, prop, varargin) % ax, lims, xrange, fmt)

            x0 = min(obj.xSurf);
            x1 = max(obj.xSurf);

            defaultAx = gca;
            defaultLims = 'auto';
            defaultLineWidth = 1.5;
            defaultFmtString = '';
            defaultXRange = [x0-1 x1+1];

            p = inputParser;

%             addRequired(p, 'prop');
            addParameter(p, 'ax', defaultAx);
            addParameter(p, 'lims', defaultLims);
            addParameter(p, 'xrange', defaultXRange);
            addParameter(p, 'fmt', defaultFmtString);
            addParameter(p, 'LineWidth', defaultLineWidth);

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

        function plt = plot_y_profile(obj, x, prop, ax)
            if nargin < 4 || isempty(ax)
                ax = gca;
                disp('Creating axes')
            end

            [q, i] = BLprof(obj,x,prop);
            size(q);
            plot(ax, q, obj.yBL(i,:))
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
            set(gca,'FontSize',12)
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

            p = inputParser;

%             addRequired(p, 'prop');
            addParameter(p, 'ax', defaultAx);
            addParameter(p, 'ploteq', false)
            addParameter(p, 'xrange', defaultXRange);
            addParameter(p, 'fmt', defaultFmtString);
            addParameter(p, 'LineWidth', defaultLineWidth);
            addParameter(p, 'LineColor', '');

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
            xlabel('H_{incomp}')
            ylabel('Pr')
            set(gca,'FontSize',12)
            disp('')
            C = colororder;
            
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
            for i=1:size(obj.yBL,1)
                Prprof = Prnow(i,1:inds(i));
                Ue = Unow(i,inds(i));
                roe = ronow(i, inds(i));
                ys = obj.yBL(i,1:inds(i));
                value(i) = trapz(ys, Prprof)/(roe*Ue^3);
            end
            value = smooth(value,obj.nSmooth);
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
            for i=1:length(inds)
                value(i) = ronow(i,inds(i)).*Unow(i,inds(i)).*(thnow(i))./munow(i,inds(i));
            end

        end

        function value = get.Ue(obj)
            inds = obj.BLedgeInd;
            Unow = obj.U;
            Misen = obj.Msurf;
            Tisen = obj.T0in./(1+((obj.gas.gam-1)/2)*Misen.^2);
            value = Misen.*sqrt(obj.gas.cp*(obj.gas.gam-1)*Tisen);                                               
%             for i=1:length(inds)
%                 value(i) = Unow(i,inds(i));
%             end
        end

        function value = get.cf(obj)
            inds = obj.BLedgeInd;
            ronow = obj.oGridProp('ro');
            Unow = obj.U;
            for i=1:length(inds)
                roe(i) = ronow(i,inds(i));
                Ue(i) = Unow(i,inds(i));
            end
                value = obj.tau_w'./(0.5*roe.*Ue.*Ue);
        end


        
        function value = get.ctau(obj)
            value = obj.get_ctau;
        end

        function value = get_ctau(obj)
            inds = obj.BLedgeInd;
            Uenow = obj.Ue;
            ronow = obj.oGridProp('ro');
            for i=1:length(inds)
                roe(i) = ronow(i,inds(i));
            end
            
            Prnow = obj.oGridProp('Pr');
            dUdynow = obj.dUdy;
            tau = Prnow./dUdynow;
            value = tau./(roe'.*Uenow.^2);
        end

        function value = get.ctau_max(obj)
            ctau = obj.ctau;
            value = max(ctau,[],2);
        end
                
        function value = y_ctau_max(obj)
            ctau = obj.ctau;
            [~, inds] = max(ctau,[],2);
            for i = 1:length(inds)
                value(i) = obj.yBL(i,inds(i));
            end
        end

        function value = j_ctau_max(obj)
            ctau = obj.ctau;
            [~, value] = max(ctau,[],2);
        end

                

        function value = get.pdyn(obj)
            inds = obj.BLedgeInd;
            ronow = obj.oGridProp('ro');
            Unow = obj.U;
            for i=1:length(inds)
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

            value = sqrt(obj.tau_w./ro_w);

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


        
        function blContour(obj, prop, ax, lims, fmt)
            if nargin < 3 || isempty(ax)
                ax = gca;
            end
            q = obj.oGridProp(prop);

            if nargin > 4 && ~isempty(fmt)
                pcolor(ax,obj.xO,obj.yO,q,fmt);
            else
                pcolor(ax,obj.xO,obj.yO,q);
            end
            shading interp
            axis equal
            if nargin > 3 && ~isempty(lims)
                caxis(lims);
            end

            cb = colorbar(ax,"southoutside");
            cb.Label.String = string(prop);

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

        
    end
end
