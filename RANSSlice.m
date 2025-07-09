classdef RANSSlice < aveSlice
    % MEANSLICE Contains the 2D (spanwise averaged) mean flow

        
    properties
        %blk;
        turb_model;
        solver;
        StR_store;            % Strain rate
        mut_store;            % Eddy viscosity
        mut_ratio_store;      % Eddy viscosity ratio
        Pr_store;
        tau_store;
%         StR;
%         mut;
        k;                    % TKE
        omega;                % Specific dissipation rate
        sp;                   % Spalart variable
        Reth;
        gamma;
        trans;          % Transition model on/off
        mod;            % beta1 modification on/off
        beta1_fac;
        metadata;
    end

    properties (Dependent = true, Hidden = true)
        Pr;             % Turbulence production
        Pr_dist;
        tau_Re;
    end

    methods
        function obj = RANSSlice(blk, gas, bcs, solver, path, iPS)
            if nargin < 6
                iPS = false;
            end
            obj@aveSlice(blk, gas, bcs, iPS);
            disp('Constructing RANSSlice')


            if nargin > 3 && ~isempty(solver)

                obj.solver = solver;
                switch solver
                    case '3dns'

                        % Get latest mean_count
                        [~, folder] = fileparts(path);
                        if strcmp(string(folder(1:3)), "run")
                            fid = fopen(fullfile(path, 'mean_flo', 'mean_time.txt'),'r');
                        else
                            fid = fopen(fullfile(path, 'mean_time.txt'),'r');
                        end
                        while ~feof(fid) % Use lastest mean files
                            temp=fgetl(fid);
                        end
                        fclose(fid);
                        temp = str2num(temp);
                        mean_count = temp(1);

                        for ib = 1:obj.NB

                            ni = blk.blockdims(ib,1);
                            nj = blk.blockdims(ib,2);
                            nk = blk.blockdims(ib,3);
                            
                            flowpath = fullfile(path, sprintf('flow2_%d_%d', [ib mean_count]));
                            f = dir(flowpath);
                            flofile = fopen(flowpath);
                            ransfile = fopen(fullfile(path, sprintf('rans_%d', ib)));
                            
                            nstats = f.bytes/(ni*nj*8);
                            A = fread(flofile,ni*nj*nstats,'float64');
                            A = reshape(A,nstats,[])';
                            fclose(flofile);

                            obj.ro{ib} = reshape(A(:,1),ni,nj);
                            obj.u{ib} = reshape(A(:,2),ni,nj)./obj.ro{ib};
                            obj.v{ib} = reshape(A(:,3),ni,nj)./obj.ro{ib};
                            obj.w{ib} = reshape(A(:,4),ni,nj)./obj.ro{ib};
                            obj.Et{ib} = reshape(A(:,5),ni,nj);

                            A = fread(ransfile,ni*nj*nk,'float64');
                            A = reshape(A,ni,nj,nk);
                            fclose(ransfile);
                            
                            obj.mut_ratio_store{ib} = mean(A,3);
                            

                        end
                end

            end
        end



        function contours = kPlot(obj,prop,ax,lims,label)
            disp('here')
            
            if nargin < 3 || isempty(ax)
                ax = gca;
            end
            q = obj.(prop);
            hold on
            for i=1:obj.NB
                contours{i} = pcolor(ax, obj.blk.x{i}, obj.blk.y{i}, q{i});
            end
            shading('interp')
%             pbaspect([6 2 1])
%             axis([0.3 0.9 0 0.2])
            %axis equal
            %axis off
            cb = colorbar('southoutside');
            if prop == 'M'
                cb.Label.String = 'M';
                cb.Label.Interpreter = 'latex';
            end
            if nargin > 3 && ~isempty(lims)
                caxis(lims)
            end
            if nargin > 4 && ~isempty(label)
                cb.Label.String = label;
                cb.Label.Interpreter = 'latex';
            end
            axis equal
            if ~isempty(obj.blk.viewarea)
                aspect = [(obj.blk.viewarea(2)-obj.blk.viewarea(1)) (obj.blk.viewarea(4)-obj.blk.viewarea(3)) 1];
                pbaspect(aspect)
                axis(obj.blk.viewarea);
            end

        end

        function value = get.tau_Re(obj)

            value = cell(1,obj.NB);
            St = obj.St_an;
            deltaij(1,1,:,:) = [1 0 0; 0 1 0; 0 0 1];

            if isempty(obj.tau_store)
                for ib=1:obj.NB
                    value{ib} = 2*obj.mut_store{ib}.*St{ib};
                    if ~strcmp(obj.turb_model, "sa")
                        value{ib} = value{ib} - (2/3)*(obj.ro{ib}.*obj.k{ib}).*deltaij;
                    end
                end
            else
                value = obj,tau_store;
            end


        end

        function value = get.Pr(obj)
            value = cell(1,obj.NB);
            if ~isempty(obj.Pr_store)
                disp('Using stored Pr')
                value = obj.Pr_store;
            elseif ~isempty(obj.mut_store)
                disp('Calculating Pr: using stored mut')
                deltaij(1,1,:,:) = [1 0 0; 0 1 0; 0 0 1];
                St_an_now = obj.St_an;
                vgrad = obj.duidxj;
                if strcmp(obj.turb_model, "sa")
                    for nb = 1:obj.NB
                        tauij = 2*obj.mut_store{nb}.*St_an_now{nb};
                        value{nb} = sum(sum(vgrad{nb}.*tauij, 4), 3);
                    end
                else
                    for nb = 1:obj.NB
                        tauij = 2*obj.mut_store{nb}.*St_an_now{nb} - (2/3)*(obj.ro{nb}.*obj.k{nb}).*deltaij;
                        value{nb} = sum(sum(vgrad{nb}.*tauij, 4), 3);
                    end
                end
            elseif ~isempty(obj.mut_ratio_store)
                disp('Calculating Pr: using stored mut ratio')
                deltaij(1,1,:,:) = [1 0 0; 0 1 0; 0 0 1];
                St_an_now = obj.St_an;
                vgrad = obj.duidxj;
                munow = obj.mu;
                for nb = 1:obj.NB
                    tauij = 2*munow{nb}.*obj.mut_ratio_store{nb}.*St_an_now{nb} - (2/3)*(obj.ro{nb}.*obj.k{nb}).*deltaij;
                    value{nb} = sum(sum(vgrad{nb}.*tauij, 4), 3);
                end
            else
                disp('Calculating Pr with k-om SST formulation')
                value = obj.mut_koSST;
            end

        end

        function value = get_StR(obj)
            if isempty(obj.StR_store)
                value = cell(1,obj.NB);
                for nb =1:obj.NB
                    disp('Calculating strain rate magnitude')
                    value{nb} = strain_rate_magnitude(obj.blk.x{nb}, obj.blk.y{nb}, obj.u{nb}, obj.v{nb});
                end
            else
                disp('Using stored strain rate magnitude')
                value = obj.StR_store;
            end
        end

        function value = mut_koSST(obj)
            value = cell(1,obj.NB);
            StRnow = obj.StR;
            for ib = 1:obj.NB
                arg2 = max(2*sqrt(obj.k{ib})./(0.09*obj.omega{ib}.*obj.wallDist{ib}), ...
                    500*obj.nu{ib}./(obj.omega{ib}.*obj.wallDist{ib}.^2));
                F2 = tanh(arg2.^2);
                mutnow = 0.31*obj.ro{ib}.*obj.k{ib}./max(0.31*obj.omega{ib}, StRnow{ib}.*F2);
                value{ib} = mutnow.*StRnow{ib}.^2;
            end
        end

        function value = Pr_koSST(obj)
            value = cell(1,obj.NB);
        end

        function obj.set_mut(obj, value)
            obj.mut_store = value;
        end

        function value = get_mut(obj)
            if ~isempty(obj.mut_store)
                value = obj.mut_store;
            elseif ~isempty(obj.mut_ratio_store)
                munow = obj.mu;
                for ib = 1:obj.NB
                    value{ib} = munow{ib}.*obj.mut_ratio_store{ib};
                end
            elseif strcmp(obj.turb_model , 'ko')
                value = obj.mut_koSST;
            else
                
            end
        end

        function value = get_mut_ratio(obj)
            if ~isempty(obj.mut_ratio_store)
                value = obj.mut_ratio_store;
            else
                mutnow = obj.mut;
                munow = obj.mu;
                for ib = 1:obj.NB
                    value{ib} = mutnow{ib}./munow{ib};
                end
            end
        end


        function value = mut_opt_cleaned(obj,smoothwindow)

            if nargin < 2
            end

            value = obj.mut;
            
        end
            
        
    end
end