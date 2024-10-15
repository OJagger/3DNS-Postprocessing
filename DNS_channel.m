classdef DNS_channel < DNS_case
    %DNS_CHANNEL Subclass of DNS_case contaning methods and properties
    ... specific to channel flow cases

    properties
        nbi;
        nbj;
    end

    properties (Dependent = true)
        Re_theta_in
    end

    methods
        function obj = DNS_channel(casename,run)
            %DNS_CHANNEL Construct an instance of this class

            args.topology = 3;
            if nargin > 0
                if nargin < 2
                    run = [];
                end
            else
                casename = [];
                run = [];
            end
            obj@DNS_case(casename,run,args);

            if nargin > 0 && ~isempty(casename)
                obj.compute_blk_metadata;
            end
        end

        function compute_blk_metadata(obj)
            obj.blk.viewarea = [inf -inf inf -inf];
            for ib = 1:obj.NB
                obj.blk.viewarea(1) = min(obj.blk.viewarea(1), min(obj.blk.x{ib},[],'all'));
                obj.blk.viewarea(2) = max(obj.blk.viewarea(2), max(obj.blk.x{ib},[],'all'));
                obj.blk.viewarea(3) = min(obj.blk.viewarea(3), min(obj.blk.y{ib},[],'all'));
                obj.blk.viewarea(4) = max(obj.blk.viewarea(4), max(obj.blk.y{ib},[],'all'));

                obj.blk.walldist{ib} = obj.blk.y{1};
            end
            obj.blk.aspect = [(obj.blk.viewarea(2)-obj.blk.viewarea(1)) ...
                (obj.blk.viewarea(4)-obj.blk.viewarea(3)) 1];

            obj.nbi = 1;
            while obj.blk.next_block{obj.nbi}.ip ~= 0
                obj.nbi = obj.nbi+1;
            end
            obj.nbj = obj.NB/obj.nbi;
        end

        function split_domain(obj, newCase)
            newFlow = volFlow();
            newFlow.flowpath = newCase.casepath;
            if isempty(obj.instFlow)
                obj.readInstFlow;
            end
            ims = ones(newCase.nbi, newCase.nbj);
            ips = ones(newCase.nbi, newCase.nbj);
            jms = ones(newCase.nbi, newCase.nbj);
            jps = ones(newCase.nbi, newCase.nbj);

            X = obj.concat_prop(obj.blk.x);
            Y = obj.concat_prop(obj.blk.y);

            for j = 1:newCase.nbj
                ib = 1+(j-1)*newCase.nbi;
                xc = newCase.blk.x{ib}(end,end);
                yc = newCase.blk.y{ib}(end,end);
                dist = sqrt((X-xc).^2+(Y-yc).^2);
                [~, ind] = min(dist,[],'all');
                [ic, jc] = ind2sub(size(X), ind);
                jps(:,j) = jc;
            end

            for i = 1:newCase.nbi
                ib = i;
                xc = newCase.blk.x{ib}(end,end);
                yc = newCase.blk.y{ib}(end,end);
                dist = sqrt((X-xc).^2+(Y-yc).^2);
                [~, ind] = min(dist,[],'all');
                [ic, jc] = ind2sub(size(X), ind);
                ips(i,:) = ic;
            end
            jms(:,2:end) = jps(:,1:end-1);
            ims(2:end,:) = ips(1:end-1,:);

            fkc = linspace(0,1,obj.blk.nk);
            fk = linspace(0,1,newCase.blk.nk);
            props = {'ro','u','v','w','Et'};
            for ip = 1:length(props)
                prop = props{ip};
                fprintf('Interpolating %s\n', prop);
                propnow = obj.concat_prop(obj.instFlow.(prop));

                for j = 1:newCase.nbj
                    for i = 1:newCase.nbi
                        ib = i+(j-1)*newCase.nbi;
                        fprintf('Block %d\n', ib);
                        Xn = X(ims(i,j):ips(i,j),jms(i,j):jps(i,j));
                        Yn = Y(ims(i,j):ips(i,j),jms(i,j):jps(i,j));
                        Vn = propnow(ims(i,j):ips(i,j),jms(i,j):jps(i,j),:);
                        [fic, fjc] = newFlow.get_spacing(Xn,Yn);
                        [fi,fj] = newFlow.get_spacing(newCase.blk.x{ib}, newCase.blk.y{ib});
                        [Jc,Ic,Kc] = meshgrid(fjc,fic,fkc);
                        [J,I,K] = meshgrid(fj,fi,fk);
%                         newFlow.(prop){ib} = interp3(Jc,Ic,Kc,Vn,J,I,K);
                        data = interp3(Jc,Ic,Kc,Vn,J,I,K);
                        save(fullfile(newCase.casepath,sprintf('block_%d_%s.mat',[ib, prop])), 'data', '-v7.3');
                        clear data
                    end
                end

                clear propnow
            end

            for i=1:newCase.NB
                flow = volFlowBlock();
                flow.ib = i;
                for ip = 1:length(props)
                    fprintf('Reading %s, block %d\n',props{ip}, i)
                    prop = props{ip};
                    block = load(fullfile(newCase.casepath,sprintf('block_%d_%s.mat',[i, prop])));
                    flow.(prop) = block.data;
                    clear data
                end

                flow.writeFlow(newCase.casepath, newCase.casetype);

            end

            newFlow.blk = newCase.blk;
            newFlow.gas = newCase.gas;
            newFlow.NB = newCase.NB;
            newFlow.casetype = obj.casetype;
            newCase.instFlow = newFlow;
        end

        function assemble_flow(obj)
            props = {'ro','u','v','w','Et'};
            
            for i=1:obj.NB
                fprintf('Assembling block %d\n', i)
                flow = volFlowBlock();
                for ip = 1:length(props)
                    prop = props{ip};
                    fprintf('Reading %s\n', prop);
                    block = load(fullfile(obj.casepath,sprintf('block_%d_%s.mat',[i, prop])));
                    flow.(prop) = block.data;
                    clear data
                end
                flow.ib = i;
   
                fprintf('Writing flow\n')
                flow.writeFlow(obj.casepath, obj.casetype);
                clear flow
            end
        end

        function arr = concat_prop(obj, prop)
            arr = [];
            for j = 1:obj.nbj
                row = [];
                for i=1:obj.nbi
                    ib = i+(j-1)*obj.nbi;
                    row = [row; prop{ib}];
                    if i < obj.nbi
                        row = row(1:end-1,:,:);
                    end
                end
                arr = [arr row];
                if j < obj.nbj
                    arr = arr(:,1:end-1,:);
                end
            end
        end

        function newCase = instantiate(obj)
            newCase = DNS_channel;
        end

        function data = get_inlet_bc_data(obj, x)

            if nargin < 2 || isempty(x)
                is = 1;
                x = 0;
                xturb = 0.25;
            else
                is = obj.meanFlow.x2ind(x);
                xturb = x;
            end
            
            data.xsample = x;
            del = obj.meanFlow.x2prop(x, 'delta99');
            y = obj.blk.y{1}(is,:);
            data.x = obj.blk.x{1}(is,:);
            data.y = y;
            data.y(end) = obj.blk.y{1}(1,end);
            data.kprof = obj.meanFlow.BLprof(xturb,'k');
            data.mut_opt = obj.meanFlow.BLprof(xturb,'mut_opt');
            data.mut_ratio = obj.meanFlow.BLprof(xturb,'mut_opt_ratio');
            data.omprof = obj.meanFlow.BLprof(xturb,'omega_opt');
            data.uprof = obj.meanFlow.BLprof(x,'u');
            data.vprof = obj.meanFlow.BLprof(x,'v');
            Vin = sqrt(data.uprof.^2+data.vprof.^2);
            data.Vin = Vin;%blend2freestream(y, Vin, del, obj.bcs.vin);
            data.Tu = sqrt(2*data.kprof/3)./data.Vin;
            [data.Toin, ~] = obj.meanFlow.BLprof(x,'T0');
            data.Tin = data.Toin - data.Vin.^2/(2*obj.gas.cp);
            Minf = M_VT0(obj.bcs.vin,obj.bcs.Toin,obj.gas.gam,obj.gas.cp);
            Mprof = M_VT0(data.Vin,data.Toin,obj.gas.gam,obj.gas.cp);
            data.Mprof = Mprof*Minf/Mprof(end);

            data.uprof = blend2freestream(y, data.uprof, del, obj.bcs.vin);
            data.vprof = blend2freestream(y, data.vprof, del, 0);
%             data.omprof = blend2freestream(y, data.omprof, del, 100);

            ps = obj.bcs.Poin*p_p0(Minf, obj.gas.gam);
            data.pprof = ps*ones(size(data.x));
            rgas = obj.gas.cp * (1- 1/obj.gas.gam);
            data.roprof = data.pprof ./ (rgas * data.Tin);
            data.Poin = ps./p_p0(data.Mprof,obj.gas.gam);
%             uprof = obj.meanFlow.BLprof(0.25,'u');
%             vprof = obj.meanFlow.BLprof(0.25,'v');
            data.aprof = atand(data.vprof./data.uprof);
%             [Poin, ~] = obj.inletProf(obj.meanFlow,'p0');
%             [Mprof, ~] = obj.inletProf(obj.meanFlow,'M');
            data.pprof = data.Poin.*p_p0(data.Mprof,obj.gas.gam);

            function profnow = blend2freestream(y, prof, delta, freeval)

                fac = 0.5*(1+tanh(y/delta - 2));
                profnow = (1-fac).*prof + fac.*freeval;

            end

        end

        function s = write_hydra_inlet_bc(obj, x, iwrite)

            if nargin < 2 || isempty(x)
                is = 1;
                x = 0;
            end
            if ~isfloat(x)
                data  = x;
            else
                data = obj.get_inlet_bc_data(x);
            end

            if nargin < 3 || isempty(iwrite)
                iwrite = true;
            end

            loc = '';
            if data.xsample > 0
                loc = strrep(sprintf('_%4.2f', data.xsample),'.','-');
            end
            s = ['hydra_inlet_bc_data' loc '.xml'];
            s = fullfile(obj.casepath, s);

%             kprof = obj.meanFlow.BLprof(0.25,'k');
%             omprof = obj.meanFlow.BLprof(0.25,'omega_opt');
%             [uprof, y] = obj.inletProf(obj.meanFlow,'u');
%             [vprof, ~] = obj.inletProf(obj.meanFlow,'v');
%             Vin = sqrt(uprof.^2+vprof.^2);
%             [Toin, ~] = obj.inletProf(obj.meanFlow,'T0');
%             Tin = Toin - Vin.^2/(2*obj.gas.cp);
%             Minf = M_VT0(obj.bcs.vin,obj.bcs.Toin,obj.gas.gam,obj.gas.cp);
%             Mprof = M_VT0(Vin,Toin,obj.gas.gam,obj.gas.cp);
%             Mprof = Mprof*Minf/Mprof(end);
%             ps = obj.bcs.Poin*p_p0(Minf, obj.gas.gam);
%             Poin = ps./p_p0(Mprof,obj.gas.gam);
% %             uprof = obj.meanFlow.BLprof(0.25,'u');
% %             vprof = obj.meanFlow.BLprof(0.25,'v');
%             aprof = atand(vprof./uprof);
% %             [Poin, ~] = obj.inletProf(obj.meanFlow,'p0');
% %             [Mprof, ~] = obj.inletProf(obj.meanFlow,'M');
%             pprof = Poin.*p_p0(Mprof,obj.gas.gam);
            

%             [vel_prof, po_prof, To_prof, T_prof] = blasius_bl(obj.bcs.Toin, obj.bcs.vin, obj.bcs.theta, obj.blk.y{1}(1,:), obj.gas);
%             Toin = To_prof*obj.bcs.Toin;
%             Tin = obj.bcs.Toin - obj.bcs.vin^2/(2*obj.gas.cp);
%             Poin = po_prof*obj.bcs.Poin;
%             rgas = obj.gas.cp*(obj.gas.gam-1)/obj.gas.gam;
%             Mprof = obj.bcs.vin*vel_prof./sqrt(obj.gas.gam*rgas*Tin*T_prof);

            props2write = ["y" "Poin" "Toin" "aprof" "Mprof" "kprof" "omprof"];
            prop_names = ["y" "ptotal" "ttotal" "apitch" "mach" "turbke" "turbomga"];

            for np=1:length(props2write)
                prop = props2write(np);
                prof = data.(prop);
                prof(isnan(prof))=0;
                A(:,np) = prof;
            end

            if iwrite
                f = fopen(s,'w');
    
                fprintf(f, '<?xml version="1.0" encoding="UTF-8"?>\n');
                fprintf(f, '<hydra>\n');
                fprintf(f, '\t<bc>\n');
                fprintf(f, '\t\t<table rank="1" ncols="%d" units="si" form="unstructured" interp="spline">\n',size(A,2)-1);
                fprintf(f, '\t\t\t<npoints>%d</npoints>\n', size(A,1));
                fprintf(f, '\t\t\t<names>\n');
                header = repmat('%s ',1,size(A,2));
                header = ['\t\t\t\t' header(1:end-1) '\n'];
                fprintf(f, header, prop_names);
                fprintf(f, '\t\t\t</names>\n');
    
                fprintf(f, '\t\t\t<data>\n');
                fmt = repmat('%8.6f ',1,size(A,2));
                fmt = ['\t\t\t\t' fmt(1:end-1) '\n'];
                fprintf(f, fmt, reshape(A',1,[]));
                fprintf(f, '\t\t\t</data>\n');
    
                fprintf(f, '\t\t</table>\n');
                fprintf(f, '\t</bc>\n');

                fprintf(f, '</hydra>\n');
    
    %             for j=1:length(kprof)
    % %                 data = fprintf(f,'\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\n',...
    % %                     [y(j) Toin(j) Poin(j) Mprof(j) aprof(j) kprof(j) omprof(j)]);
    % %                 data = fprintf(f,'\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\n',...
    % %                     [y(j) Poin(j) Toin(j) Mprof(j) aprof(j) 0.0 1000.0]);
    %                 data = fprintf(f,'\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\n',...
    %                     [y(j) Toin(j) Poin(j) aprof(j) kprof(j) omprof(j)]);
    %             end
                fclose(f);
            end
        end

        function s = write_fluent_inlet_bc(obj, x, iwrite)

            if nargin < 2 || isempty(x)
                is = 1;
                x = 0;
                xturb = 0.25;
            elseif ~isfloat(x)
                data  = x;
            else
                data = obj.get_inlet_bc_data(x);
            end

            if nargin < 3 || isempty(iwrite)
                iwrite = true;
            end

            loc = '';
            if data.xsample > 0
                loc = strrep(sprintf('_%4.2f', data.xsample),'.','-');
            end
            s = ['fluent_inlet_bc_data' loc '.prof'];
            s = fullfile(obj.casepath, s);
    
            if iwrite
                f = fopen(s,'w');
    
                props2write = ["x" "y" "uprof" "vprof" "Vin" "Tin" "pprof" "Toin" "Poin" "roprof" "aprof" "kprof" "omprof" "mut_opt" "mut_ratio" "Tu" "Mprof"];
                prop_names = ["x" "y" "u" "v" "vmag" "t" "p" "tstag" "pstag" "ro" "a" "k" "omega" "mut_opt" "mut_ratio" "tu" "mach"];
    
                fprintf(f, '((inlet-bl line %d)', length(data.x));
                for np = 1:length(prop_names)
                    prop = props2write(np);
                    prof = data.(prop);
                    prof(isnan(prof))=0;
    
                    fprintf(f, '\n(%s', prop_names(np));
                    fprintf(f, "\n%f %f %f %f %f %f %f %f %f %f %f %f", prof);
                    fprintf(f, "\n)");
                end
                fprintf(f, ')');
                fclose(f);
            end

        end

        function [s, ps] = write_fluent_freestream_bc(obj, iwrite)

            if nargin < 2 || isempty(iwrite)
                iwrite = true;
            end

            data.x = obj.blk.x{1}(:,end);
            data.y = obj.blk.y{1}(:,end);
            data.pprof = obj.meanFlow.p{1}(:,end);
            ps = data.pprof(end);

            s = fullfile(obj.casepath, 'fluent_freestream_bc_data.prof');

            if iwrite
                f = fopen(s,'w');

                props2write = ["x" "y" "pprof"];
                prop_names = ["x" "y" "p"];

                fprintf(f, '((freestream line %d)', length(data.x));
                for np = 1:length(prop_names)
                    prop = props2write(np);
                    prof = data.(prop);
                    prof(isnan(prof))=0;
    
                    fprintf(f, '\n(%s', prop_names(np));
                    fprintf(f, "\n%f %f %f %f %f %f %f %f %f %f %f %f", prof);
                    fprintf(f, "\n)");
                end
                fprintf(f, ')');
                fclose(f);
            end

        end

        function [s, ps] = write_hydra_freestream_bc(obj, iwrite)

            if nargin < 2 || isempty(iwrite)
                iwrite = true;
            end

            data.x = obj.blk.x{1}(:,end);
            data.y = obj.blk.y{1}(:,end);
            data.pprof = obj.meanFlow.p{1}(:,end);
            ps = data.pprof(end);

            s = 'hydra_freestream_bc_data.xml';
            s = fullfile(obj.casepath, s);

            props2write = ["x" "pprof"];
            prop_names = ["x" "pstatic"];

            for np=1:length(props2write)
                prop = props2write(np);
                prof = data.(prop);
                prof(isnan(prof))=0;
                A(:,np) = prof;
            end
    
            if iwrite
                f = fopen(s,'w');
    
                fprintf(f, '<?xml version="1.0" encoding="UTF-8"?>\n');
                fprintf(f, '<hydra>\n');
                fprintf(f, '\t<bc>\n');
                fprintf(f, '\t\t<table rank="1" ncols="%d" units="si" form="unstructured" interp="spline">\n',size(A,2)-1);
                fprintf(f, '\t\t\t<npoints>%d</npoints>\n', size(A,1));
                fprintf(f, '\t\t\t<names>\n');
                header = repmat('%s ',1,size(A,2));
                header = ['\t\t\t\t' header(1:end-1) '\n'];
                fprintf(f, header, prop_names);
                fprintf(f, '\t\t\t</names>\n');
    
                fprintf(f, '\t\t\t<data>\n');
                fmt = repmat('%8.6f ',1,size(A,2));
                fmt = ['\t\t\t\t' fmt(1:end-1) '\n'];
                fprintf(f, fmt, reshape(A',1,[]));
                fprintf(f, '\t\t\t</data>\n');
    
                fprintf(f, '\t\t</table>\n');
                fprintf(f, '\t</bc>\n');

                fprintf(f, '</hydra>\n');
    
                fclose(f);
            end

        end

        function write_hydra_split_inlet_bc(obj)
            kprof = obj.meanFlow.BLprof(0.25,'k');
            omprof = obj.meanFlow.BLprof(0.25,'omega_opt');
            [uprof, y] = obj.inletProf(obj.meanFlow,'u');
            [vprof, ~] = obj.inletProf(obj.meanFlow,'v');
            Vin = sqrt(uprof.^2+vprof.^2);
            [Toin, ~] = obj.inletProf(obj.meanFlow,'T0');
            Tin = Toin - Vin.^2/(2*obj.gas.cp);
            Minf = M_VT0(obj.bcs.vin,obj.bcs.Toin,obj.gas.gam,obj.gas.cp);
            Mprof = M_VT0(Vin,Toin,obj.gas.gam,obj.gas.cp);
            Mprof = Mprof*Minf/Mprof(end);
            ps = obj.bcs.Poin*p_p0(Minf, obj.gas.gam);
            Poin = ps./p_p0(Mprof,obj.gas.gam);
%             uprof = obj.meanFlow.BLprof(0.25,'u');
%             vprof = obj.meanFlow.BLprof(0.25,'v');
            aprof = atand(vprof./uprof);
%             [Poin, ~] = obj.inletProf(obj.meanFlow,'p0');
%             [Mprof, ~] = obj.inletProf(obj.meanFlow,'M');
            pprof = Poin.*p_p0(Mprof,obj.gas.gam);
            

%             [vel_prof, po_prof, To_prof, T_prof] = blasius_bl(obj.bcs.Toin, obj.bcs.vin, obj.bcs.theta, obj.blk.y{1}(1,:), obj.gas);
%             Toin = To_prof*obj.bcs.Toin;
%             Tin = obj.bcs.Toin - obj.bcs.vin^2/(2*obj.gas.cp);
%             Poin = po_prof*obj.bcs.Poin;
%             rgas = obj.gas.cp*(obj.gas.gam-1)/obj.gas.gam;
%             Mprof = obj.bcs.vin*vel_prof./sqrt(obj.gas.gam*rgas*Tin*T_prof);
            
            [~, jsonic] = min(abs(Mprof-1));
            
            f = fopen(fullfile(obj.casepath, 'hydra_inlet_bc_sub.txt'),'w');
            for j=1:jsonic
                data = fprintf(f,'\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\n',...
                    [y(j) Poin(j) Toin(j) aprof(j) kprof(j) omprof(j)]);
            end
            fclose(f);
            
            f = fopen(fullfile(obj.casepath, 'hydra_inlet_bc_sup.txt'),'w');
            for j=jsonic:length(kprof)
                data = fprintf(f,'\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\t%8.6f\n',...
                    [y(j) Poin(j) Toin(j) Mprof(j) aprof(j) kprof(j) omprof(j)]);
            end
            fclose(f);
        end

        function boundaries = getBoundaries(obj)


            % BCs:
            % 3: Wall
            % 4: Pressure inlet
            % 5: Pressure outlet
            % 12: Periodic
            % 8: Periodic shadow

            boundaries = {};
            b.label = "Inlet";
            blks = [];
            for j = 1:obj.nbj
                ib = 1+(j-1)*obj.nbi;
                blks = [blks ib];
            end
            b.blocks = blks;
            b.patches(1:length(blks)) = 1;
            b.type = 4;
            boundaries{end+1} = b;
            
            b.label = "Wall";
            b.blocks = 1:obj.nbi;
            b.patches(1:obj.nbi) = 3;
            b.type = 3;
            boundaries{end+1} = b;

            b.Label = "Freestream";
            blks = (obj.NB-obj.nbi)+1:obj.NB;
            b.patches(1:length(blks)) = 4;
            b.type = 5;
            boundaries{end+1} = b;
            
%             b.label = "Pre-shock";
%             xmid = obj.blk.x{obj.nbi}(end,1)/2;
%             blks = [];
%             psblks = [];
%             for ib=(obj.NB-obj.nbi)+1:obj.NB
%                 if obj.blk.x{ib}(floor(obj.blk.blockdims(ib,1)/2),end) < xmid
%                     blks = [blks ib];
%                 else
%                     psblks = [psblks ib];
%                 end
%             end
%             b.blocks = blks;
%             b.patches(1:length(blks)) = 4;
%             b.type = 5;
%             boundaries{end+1} = b;
% 
%             b.label = "Post-shock";
%             b.blocks = psblks;
%             b.patches(1:length(psblks)) = 4;
%             b.type = 5;
%             boundaries{end+1} = b;
            
            b.label = "Outlet";
            blks = [];
            for j = 1:obj.nbj
                ib = j*obj.nbi;
                blks = [blks ib];
            end
            b.blocks = blks;
            b.patches(1:length(blks)) = 2;
            b.type = 5;
            boundaries{end+1} = b;
        end

        function write_3DNS_inlet_bc(obj, x, path, blk, iwrite)

            if nargin < 2 || isempty(x)
                data = obj.get_inlet_bc_data;
            elseif ~isfloat(x)
                data  = x;
            else
                data = obj.get_inlet_bc_data(x);
            end

            if nargin < 3 || isempty(path)
                path = obj.casepath;
            end

            if nargin < 4 || isempty(blk)
                blk = obj.blk;
            end

            if nargin < 5 || isempty(iwrite)
                iwrite = true;
            end

            To = interp1(data.y, data.Toin/data.Toin(end), blk.y{1}(1,:));
            Po = interp1(data.y, data.Poin/data.Poin(end), blk.y{1}(1,:));
            vel = interp1(data.y, data.Vin/data.Vin(end), blk.y{1}(1,:));
            a = interp1(data.y, pi*data.aprof/180, blk.y{1}(1,:));

            if iwrite
                f = fopen(fullfile(path, 'inlet_profile.txt'), 'w');
                fprintf(f, '%4.2f\n', 0.0);

                fprintf(f, '%6.6f %8.6f %8.6f %8.6f\n', [vel; Po; To; a]);

                fclose(f);
            end

        end
        
        function write_3DNS_freestream_bc(obj, Min, xShock, Lshock, path)

            if nargin < 5
                path = obj.casepath;
            end

            gam = obj.gas.gam;
            cp = obj.gas.cp;
            rgas = cp*(gam-1)/gam;

            % Pre-shock conditions
            fM = 1+0.5*(gam-1)*Min^2;
            pin = obj.bcs.Poin*fM^(-gam/(gam-1));
            tin = obj.bcs.Toin/fM;
            roin = pin/(rgas*tin);
            vin = Min*sqrt(gam*rgas*tin);

            % Post shock conditions
            Ms = sqrt(fM/(gam*Min^2 - 0.5*(gam-1)));
            ps = pin*(1+2*gam*(Min^2-1)/(gam+1));
            ros = 0.5*roin*(gam+1)*Min^2/fM;
            Ts = ps/(ros*rgas);
            vs = Ms*sqrt(gam*rgas*Ts);

            f = fopen(fullfile(path, 'freestream.txt'), 'w');
            fprintf(f, '%d\n', 4);
            fprintf(f, '%6.4f %8.6f %8.6f\n', [0.0 vin pin])
            fprintf(f, '%6.4f %8.6f %8.6f\n', [xShock-Lshock/2 vin pin]);
            fprintf(f, '%6.4f %8.6f %8.6f\n', [xShock+Lshock/2 vs ps]);
            fprintf(f, '%6.4f %8.6f %8.6f\n', [1.0 vs ps]);

        end

        function [flow vin ps muref] = init_shock_flow(obj, Min, xShock, Lshock, theta_in, Reth_in)
            gam = obj.gas.gam;
            cp = obj.gas.cp;
            rgas = cp*(gam-1)/gam;
            flow = volFlow;
            flow.NB = 1;
            flow.gam = obj.gas.gam;
            flow.cp = obj.gas.cp;
            flow.blk = obj.blk;
            flow.gas = obj.gas;
            flow.nk = obj.solver.nk;
            flow.flowpath = obj.casepath;
            flow.casetype = 'gpu';
            [ni, nj] = size(obj.blk.x{1});

            % Pre-shock conditions
            fM = 1+0.5*(gam-1)*Min^2;
            pin = obj.bcs.Poin*fM^(-gam/(gam-1));
            tin = obj.bcs.Toin/fM;
            roin = pin/(rgas*tin);
            vin = Min*sqrt(gam*rgas*tin);
            Etin = pin/(gam-1) + 0.5*roin*vin^2;

            % Post shock conditions
            Ms = sqrt(fM/(gam*Min^2 - 0.5*(gam-1)));
            ps = pin*(1+2*gam*(Min^2-1)/(gam+1));
            ros = 0.5*roin*(gam+1)*Min^2/fM;
            Ts = ps/(ros*rgas);
            vs = Ms*sqrt(gam*rgas*Ts);
            Ets = ps/(gam-1) + 0.5*ros*vs^2;

            % Reference viscosiy
            muin = roin*vin*theta_in / Reth_in;
            muref = sutherland_mu_ref(muin, tin);

            % BL profiles
            if theta_in > 0
                [vel_prof, po_prof, To_prof, T_prof] = blasius_bl(obj.bcs.Toin, vin, theta_in, obj.blk.y{1}(1,:), obj.gas);
            else
                vel_prof = ones(1, nj);
                po_prof = ones(1, nj);
                To_prof = ones(1, nj);
                T_prof = ones(1, nj);
            end

            Et_prof = (roin./T_prof) .* ((obj.gas.cp/obj.gas.gam) * tin * T_prof + vin^2 * vel_prof.^2/2);
            Et_prof = Et_prof/Et_prof(end);


            blfn_vel = ones(ni,nj).*vel_prof;
            blfn_ro = ones(ni, nj)./T_prof;
            blfn_et = ones(ni,nj).*Et_prof;

            shfn = tanh((flow.blk.x{1}-xShock)/Lshock);


            nk = obj.solver.nk;
            flow.v{1} = zeros(ni, nj, nk);
            flow.w{1} = flow.v{1};
            flow.u{1} = 0.5*((vin + vs) - (vin-vs)*shfn).*blfn_vel;
            flow.ro{1} = 0.5*((roin + ros) - (roin-ros)*shfn).*blfn_ro;
            flow.Et{1} = 0.5*((Etin+Ets) - (Etin-Ets)*shfn).*blfn_et;


            fprintf('Inlet pressure: %7.5e\n', pin);
            fprintf('Inlet velocity: %6.2f\n', vin);            
            fprintf('Inlet temperature: %6.2fr2\n', tin);
            fprintf('Inlet density: %6.4f\n', roin);
            
            fprintf('\n')

            fprintf('Outlet pressure: %7.5e\n', ps);
            fprintf('Outlet velocity: %6.2f\n', vs);
            fprintf('Outlet Mach number: %6.4f\n', Ms);
            
        end

        function calc_muref(obj, Reth)
            
        end

        function flow = transform_flow_m1(obj, M1new)

            gam = obj.gas.gam;
            cp = obj.gas.cp;
            rgas = cp*(1-1/gam);
            rostag = obj.bcs.Poin/(rgas*obj.bcs.Toin);

            M1old = M_VT0(obj.bcs.vin, obj.bcs.Toin, obj.gas.gam, obj.gas.cp);
            fM = 1+0.5*(gam-1)*M1old^2;
            Mso = sqrt(fM/(gam*M1old^2 - 0.5*(gam-1)));
            pino = obj.bcs.Poin*fM^(-gam/(gam-1));
            pso = pino*(1+2*gam*(M1old^2-1)/(gam+1));

            fM = 1+0.5*(gam-1)*M1new^2;
            Msn = sqrt(fM/(gam*M1new^2 - 0.5*(gam-1)));
            pinn = obj.bcs.Poin*fM^(-gam/(gam-1));
            psn = pinn*(1+2*gam*(M1new^2-1)/(gam+1));

            flow = volFlow;
            flow.NB = obj.NB;
            flow.gam = obj.gas.gam;
            flow.cp = obj.gas.cp;
            flow.blk = obj.blk;
            flow.gas = obj.gas;
            flow.nk = obj.solver.nk;
            flow.flowpath = obj.casepath;
            flow.casetype = 'gpu';


            for ib = 1:obj.NB
                
%                 Mnew = Msn + (obj.instFlow.M{ib} - Mso)*(M1new - Msn)/(M1old - Mso);
                Mnew = interp1([0 Mso M1old], [0 Msn M1new], obj.instFlow.M{1}, 'linear','extrap');
                Pnew = psn + (obj.instFlow.p{ib} - pso)*(pinn - psn)/(pino - pso);
                Vnew = Vel_M(Mnew, obj.bcs.Toin, cp, gam);
                Vold = max(1e-6, obj.instFlow.vel{ib});
                Tnew = obj.bcs.Toin*T_T0(Mnew, gam);
                flow.ro{ib} = Pnew ./ (rgas * Tnew);
                flow.u{ib} = Vnew.*obj.instFlow.u{ib}./Vold;
                flow.v{ib} = Vnew.*obj.instFlow.v{ib}./Vold;
                flow.w{ib} = Vnew.*obj.instFlow.w{ib}./Vold;
                flow.Et{ib} = Pnew/(gam - 1) + 0.5*flow.ro{ib}.*Vnew.^2;

            end

        end

        function setup_freestream_pressure_dist(obj, x, u, p)

            fid = fopen(fullfile(obj.casepath, 'freestream.txt'), 'w');
            fprintf(fid,'%d\n', length(x))
            for i = 1:length(x)
                fprintf(fid, '%f %f %f\n', x(i), u(i), p(i));
            end
            fclose(fid);
        end

        function [in out] = get_BCs(obj, Min)
            gam = obj.gas.gam;
            cp = obj.gas.cp;
            rgas = cp*(gam-1)/gam;

            % Pre-shock conditions
            fM = 1+0.5*(gam-1)*Min^2;
            pin = obj.bcs.Poin*fM^(-gam/(gam-1));
            tin = obj.bcs.Toin/fM;
            roin = pin/(rgas*tin);
            vin = Min*sqrt(gam*rgas*tin);
            Etin = pin/(gam-1) + 0.5*roin*vin^2;

            in.v = vin;
            in.M = Min;
            in.T = tin;
            in.ro = roin;
            in.p = pin;

            % Post shock conditions
            Ms = sqrt(fM/(gam*Min^2 - 0.5*(gam-1)));
            ps = pin*(1+2*gam*(Min^2-1)/(gam+1));
            ros = 0.5*roin*(gam+1)*Min^2/fM;
            Ts = ps/(ros*rgas);
            vs = Ms*sqrt(gam*rgas*Ts);
            Ets = ps/(gam-1) + 0.5*ros*vs^2;

            out.v = vs;
            out.M = Ms;
            out.T = Ts;
            out.ro = ros;
            out.p = ps;

        end

        function write_mut_opt_file(obj, smoothwindow)

            if nargin < 2
                smoothwindow = 1;
            end                
           
            mto = obj.meanFlow.mut_opt_cleaned(smoothwindow);
            mu = obj.meanFlow.mu;

            for ib = 1:obj.NB
            
                f = fopen(fullfile(obj.casepath, ['mut_opt_rans_' num2str(ib)]), 'wb');
                A = reshape(mto{ib}./mu{ib},1,[]);
                fwrite(f,A,'float64');
                fclose(f);

                f = fopen(fullfile(obj.casepath, ['mut_opt2_rans_' num2str(ib)]), 'wb');
                A = zeros(1,prod(obj.blk.blockdims(ib,:)))  ;
                fwrite(f,A,'float64');
                fclose(f);
    
            end

        end

        function plotYWall(obj)

            dy = [];
            s = [];

            for i = 1:length(obj.blk.oblocks)
                ib = obj.blk.oblocks(i);
                flip = obj.blk.oblocks_flip(i);
                if obj.blk.next_patch{ib}.jm == 3
                    xwall = obj.blk.x{ib}(:,1);
                    x1 = obj.blk.x{ib}(:,2);
                    ywall = obj.blk.y{ib}(:,1);
                    y1 = obj.blk.y{ib}(:,2);
                else
                    xwall = obj.blk.x{ib}(:,end);
                    x1 = obj.blk.x{ib}(:,end-1);
                    ywall = obj.blk.y{ib}(:,end);
                    y1 = obj.blk.y{ib}(:,end-1);
                end
                
                dynow = sqrt((xwall-x1).^2 + (ywall-y1).^2);
                snow = curve_length(xwall, ywall);
                if flip
                    snow = flip(snow);
                    dynow = flip(dynow);
                end

                dynow = reshape(dynow,1,[]);

                dy = [dy dynow];
                if isempty(s)
                    s = snow;
                else
                    s = [s snow+s(end)];
                end
                    
            end

            plot(s, dy)

        end

    end
end