classdef inletProfile < handle

    properties
        blk;
        gas;
        bcs;
        casepath;
        vel_prof;
        po_prof;
        To_prof;
        alpha;
        nrelax;
        urf;
        y;
        nj;
    end

    properties (Dependent = true, Hidden = true)
        M;
        p;
        T;
        u;
        v;
        ro;
        vel;
        U;
        p0;
        T0;
        del995;
        jEdge;
        delStar;
        theta;
        delStar_k;
        theta_k;
        H;
        H_k;
    end
        
    methods
        function obj = inletProfile(blk, gas, bcs, casepath)

            if nargin > 0

                obj.blk = blk;
                obj.gas = gas;
                obj.bcs = bcs;
                
                if nargin > 3
                    obj.casepath = casepath;
                

                    f = fopen(fullfile(casepath, 'inlet_profile.txt'), 'r');
                    temp = str2num(char(split(fgetl(f))));
                    obj.nrelax = temp(1);
                    obj.urf = temp(2);
                    
                    obj.y = blk.y{blk.inlet_blocks{1}}(1,:);
                    obj.nj = length(obj.y);
    
                    data = fscanf(f, '%f %f %f %f\n', [4 Inf]);
                    obj.vel_prof = data(1,:);
                    obj.po_prof = data(2,:);
                    obj.To_prof = data(3,:);
                    obj.alpha = data(4,:);
    
                    fclose(f);
                
                end
            end
        end

        function value = get.u(obj)
            value = obj.vel.*cos(obj.alpha);
        end

        function value = get.v(obj)
            value = obj.vel.*sin(obj.alpha);
        end

        function value = get.vel(obj)
            value = obj.bcs.vin*obj.vel_prof;
        end

        function value = get.U(obj)
            value = obj.bcs.vin*obj.vel_prof;
        end

        function value = get.p0(obj)
            value = obj.bcs.Poin*obj.po_prof;
        end

        function value = get.p(obj)
            value = obj.p0.*p_p0(obj.M, obj.gas.gam);
        end

        function value = get.T0(obj)
            value = obj.bcs.Toin*obj.To_prof;
        end

        function value = get.T(obj)
            value = obj.T0.*T_T0(obj.M, obj.gas.gam);
        end

        function value = get.M(obj)
            Po = obj.bcs.Poin*obj.po_prof;
            To = obj.bcs.Toin*obj.To_prof;
            value = M_VT0(obj.vel, To, obj.gas.gam, obj.gas.cp);
        end

        function value = get.ro(obj)
            rgas = obj.gas.cp * (obj.gas.gam - 1)/obj.gas.gam;
            value = obj.p./obj.T/rgas;
        end


        function value = get.jEdge(obj)
            unow = obj.u;
            Ue = max(unow);
            jEdge = 1;
            while unow(jEdge) < Ue*0.995
                jEdge = jEdge+1;
            end
            value = jEdge;
        end

        function value = get.del995(obj)
            value = obj.y(obj.jEdge);
        end

        function value = get.theta(obj)
            
            ye = 1.1*obj.del995;
            [~, jLim] = min(abs(obj.y - ye));

            Uprof = obj.u(1:jLim);
            roprof = obj.ro(1:jLim);
            [Ue, je] = max(Uprof);
            roe = roprof(je);

            integrand = (roprof.*Uprof/(roe*Ue)).*(1-Uprof/Ue);
            value = trapz(obj.y(1:jLim), integrand);

        end

        function value = get.theta_k(obj)
            
            ye = 1.1*obj.del995;
            [~, jLim] = min(abs(obj.y - ye));

            Uprof = obj.u(1:jLim);
            [Ue, je] = max(Uprof);

            integrand = (Uprof/Ue).*(1-Uprof/Ue);
            value = trapz(obj.y(1:jLim), integrand);

        end

        function value = get.delStar(obj)
            
            ye = 1.1*obj.del995;
            [~, jLim] = min(abs(obj.y - ye));

            Uprof = obj.u(1:jLim);
            roprof = obj.ro(1:jLim);
            [Ue, je] = max(Uprof);
            roe = roprof(je);

            integrand = 1- roprof.*Uprof/(roe*Ue);
            value = trapz(obj.y(1:jLim), integrand);

        end

        function value = get.delStar_k(obj)
            
            ye = 1.1*obj.del995;
            [~, jLim] = min(abs(obj.y - ye));

            Uprof = obj.u(1:jLim);
            [Ue, je] = max(Uprof);

            integrand = 1 - Uprof/Ue;
            value = trapz(obj.y(1:jLim), integrand);

        end

        function value = get.H(obj)
            value = obj.delStar/obj.theta;
        end

        function value = get.H_k(obj)
            value = obj.delStar_k/obj.theta_k;
        end
            

        function plot(obj, prop, varargin)

            defaultAx = gca;
            defaultEdgeY = [];
            defaultLineWidth = 1.5;
            defaultFmtString = '';

            p = inputParser;

	    if ischar(prop)
                q = obj.(prop);
            else
                q = prop;
            end


%            addRequired(p, 'prop');
            addParameter(p, 'ax', defaultAx);
            addParameter(p, 'yEdge', defaultEdgeY);
            addParameter(p, 'normalise', false);
            addParameter(p, 'scale', 1);
            addParameter(p, 'fmt', defaultFmtString);
            addParameter(p, 'LineWidth', defaultLineWidth);
            addParameter(p, 'normaliseY', true);

            parse(p, varargin{:});  

            if ~isempty(p.Results.yEdge)
                yEdge = p.Results.yEdge;
            elseif p.Results.normaliseY
                yEdge = obj.del995;
            else
                yEdge = 1;
            end

            if p.Results.normalise
                scale = max(abs(q));
            else
                scale = p.Results.scale;
            end


            if isscalar(q)
                s = yline(p.Results.ax, q/yEdge, p.Results.fmt);
            else
                s = plot(p.Results.ax, q/scale, obj.y/yEdge, ...                       % corresponding y vals
                    p.Results.fmt, ...                                                                      % Set format string
                    'LineWidth', p.Results.LineWidth);
            end

        end

        function writeProf(obj, path)

            f = fopen(fullfile(path, 'inlet_profile.txt'), 'w');

            fprintf(f,'%d %8.6f\n', [obj.nrelax, obj.urf]);
            fprintf(f,'%8.6f %8.6f %8.6f %8.6f\n', [obj.vel_prof; obj.po_prof; obj.To_prof; obj.alpha]);
            fclose(f)

        end

        function prof = scale_thickness(obj, thetanow)

            prof = obj;
            prof.casepath = [];

            scale = obj.theta/thetanow;

            prof.vel_prof = interp1(obj.y/scale, obj.vel_prof, obj.y, 'pchip', 'extrap');
            prof.po_prof = interp1(obj.y/scale, obj.po_prof, obj.y, 'pchip', 'extrap');
            prof.To_prof = interp1(obj.y/scale, obj.To_prof, obj.y, 'pchip', 'extrap');
            prof.alpha = interp1(obj.y/scale, obj.alpha, obj.y, 'pchip', 'extrap');

        end

        function newprof = interpolateOntoNewCase(obj, newcase)
            
            newprof = inletProfile(newcase.blk, newcase.gas, newcase.bcs);
            newprof.y = newcase.blk.y{newcase.blk.inlet_blocks{1}}(1,:);

            newprof.nrelax = obj.nrelax;
            newprof.urf = obj.urf;
            newprof.nj = length(newprof.y);

            newprof.vel_prof = interp1(obj.y, obj.vel_prof, newprof.y, 'pchip', 'extrap');
            newprof.po_prof = interp1(obj.y, obj.po_prof, newprof.y, 'pchip', 'extrap');
            newprof.To_prof = interp1(obj.y, obj.To_prof, newprof.y, 'pchip', 'extrap');
            newprof.alpha = interp1(obj.y, obj.alpha, newprof.y, 'pchip', 'extrap');

            % newcase.inletProf = newprof;

        end

    end

end
