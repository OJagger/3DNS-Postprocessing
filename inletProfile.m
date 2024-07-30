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
        u;
        v;
        ro;
        vel;
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

        function value = get.M(obj)

            Po = obj.bcs.Toin*obj.po_prof;
            To = obj.bcs.Toin*obj.To_prof;

            value = M_VT0(obj.vel, To, obj.gas.gam, obj.gas.cp);

        end
            

        function plot(obj, prop, varargin)

            defaultAx = gca;
            defaultEdgeY = [];
            defaultLineWidth = 1.5;
            defaultFmtString = '';

            p = inputParser;

            q = obj.(prop);

            addRequired(p, 'prop', @ischar);
            addParameter(p, 'ax', defaultAx);
            addParameter(p, 'yEdge', defaultEdgeY);
            addParameter(p, 'normalise', false);
            addParameter(p, 'scale', 1);
            addParameter(p, 'fmt', defaultFmtString);
            addParameter(p, 'LineWidth', defaultLineWidth);
            addParameter(p, 'normaliseY', true);

            parse(p, prop, varargin{:});  

            if ~isempty(p.Results.yEdge)
                yEdge = p.Results.yEdge;
            else
                yEdge = 1;
            end

            if p.Results.normalise
                scale = max(abs(q));
            else
                scale = p.Results.scale;
            end


            s = plot(p.Results.ax, q/scale, obj.y/yEdge, ...                       % corresponding y vals
                p.Results.fmt, ...                                                                      % Set format string
                'LineWidth', p.Results.LineWidth);

        end

        function writeProf(obj, path)

            f = fopen(fullfile(path, 'inlet_profile.txt'), 'w');

            fprintf(f,'%d %8.6f\n', [obj.nrelax, obj.urf]);
            fprintf(f,'%8.6f %8.6f %8.6f %8.6f %8.6f\n', [obj.y; obj.vel_prof; obj.po_prof; obj.To_prof; obj.alpha]);
            fclose(f)

        end

        function interpolateOntoNewCase(obj, newcase)
            
            newprof = inletProfile(newcase.blk, newcase.gas, newcase.bcs);
            newprof.y = newcase.blk.y{newcase.blk.inlet_blocks{1}}(1,:);

            newprof.nrelax = obj.nrelax;
            newprof.urf = obj.urf;
            newprof.nj = length(newprof.y);

            newprof.vel_prof = interp1(obj.y, obj.vel_prof, newprof.y);
            newprof.po_prof = interp1(obj.y, obj.po_prof, newprof.y);
            newprof.To_prof = interp1(obj.y, obj.To_prof, newprof.y);
            newprof.alpha = interp1(obj.y, obj.alpha, newprof.y);

            newcase.inletProf = newprof;

        end

    end

end