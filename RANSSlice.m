classdef RANSSlice < aveSlice
    % MEANSLICE Contains the 2D (spanwise averaged) mean flow

        
    properties
        %blk;
        StR_store;            % Strain rate
        mut_store;            % Eddy viscosity
        Pr_store;
%         StR;
%         mut;
        k;
        omega;
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
    end

    methods
        function obj = RANSSlice(blk, gas, bcs, basedir, datafile, trans, mod, nodes, renorm)
            obj@aveSlice(blk, gas, bcs);
            disp('Constructing RANSSlice')


            if nargin > 3

                            %obj.blk = blk;
                obj.trans = trans;
                obj.mod = mod;
                transstr = ["turb","trans"];
                modstr = ["bsl","mod"];
                %data = readtable(fullfile(basedir,transstr(trans+1) + "_" + modstr(mod+1),'hf2d.txt'));
                %data = readtable(fullfile(basedir,'hf2d.txt'));
                data = readtable(fullfile(basedir,datafile));
                if any("x_coordinate" == string(data.Properties.VariableNames))
                    mode = 'fluent';
                else
                    mode = 'hydra';
                end
    
                if nargin < 7 || isempty(nodes)
                    switch mode
                        case 'fluent'
                            points.x_coordinate = data.x_coordinate;
                            points.y_coordinate = data.y_coordinate;
                            points.nodenumber = data.nodenumber;
                        case 'hydra'
                            points.x_coordinate = data.Points_0;
                            points.y_coordinate = data.Points_1;
                            points.nodenumber = 1:Psize(data,1);
                    end
                    nodes = mesh2nodes(blk,points);
                end

                if nargin < 8 || renorm == false
                    tref = 1;
                    uref = 1;
                    pref = 1;
                    roref = 1;
                    muref = 1;
                else
                    tref = 288;
                    pref = 1.013e5;
                    roref = 1.226;
                    uref = sqrt(pref/roref);
                end
                    
                
                cv = gas.cp/gas.gam;
                


                mut = cell(1,obj.NB);
                StR = cell(1,obj.NB);
                obj.wallDist = cell(1,obj.NB);
                for nb=1:obj.NB
                    fprintf('nb: %d\n', nb)
                    nib = blk.blockdims(nb,1);
                    njb = blk.blockdims(nb,2);
                    ro = zeros(nib,njb);
                    u = zeros(nib,njb);
                    v = zeros(nib,njb);
                    w = zeros(nib,njb);
                    Et = zeros(nib,njb);
                    k = zeros(nib,njb);
                    omega = zeros(nib,njb);
                    mutnow = zeros(nib,njb);
                    StRnow = zeros(nib,njb);
                    wallDistNow = zeros(nib,njb);
                    if obj.trans == 1
                        gamma = zeros(nib,njb);
                        Reth = zeros(nib,njb);
                    end

                    
                    nodesnow = nodes{nb};
                    switch mode
                        case 'fluent'
                            for i=1:nib
                                for j=1:njb
                                    node = nodesnow(i,j);
                                    rotmp = data.density(node);
                                    utmp = data.x_velocity(node);
                                    vtmp = data.y_velocity(node);
                                    Ttmp = data.temperature(node);
                                    ro(i,j) = rotmp;
                                    u(i,j) = utmp;
                                    v(i,j) = vtmp;
                                    Et(i,j) = rotmp*(cv*Ttmp + 0.5*(utmp^2+vtmp^2));
                                    k(i,j) = data.turb_kinetic_energy(node);
                                    omega(i,j) = data.specific_diss_rate(node);
                                    mutnow(i,j) = data.viscosity_turb(node);
                                    StRnow(i,j) = data.strain_rate_mag(node);
                                    if obj.trans == 1
                                        gamma(i,j) = data.intermittency(node);
                                        Reth(i,j) = data.momentum_thickness_re(node);
                                    end
                                end
                                mut{nb} = mutnow;
                                StR{nb} = StRnow;
                            end

                        case 'hydra'
                            mut=[];
                            StR = [];
                            for i=1:nib
                                for j=1:njb
                                    node = nodesnow(i,j);
                                    rotmp = data.density(node);
                                    utmp = data.relativeVelocityVector_0(node);
                                    vtmp = data.relativeVelocityVector_1(node);
                                    ptmp = data.staticPressure(node);
                                    ro(i,j) = rotmp*roref;
                                    u(i,j) = utmp*uref;
                                    v(i,j) = vtmp*uref;
                                    Ttmp = tref*ptmp/rotmp;
                                    Et(i,j) = roref*rotmp*(cv*Ttmp + 0.5*uref^2*(utmp^2+vtmp^2));
                                    k(i,j) = data.turbulentKineticEnergy(node);
                                    omega(i,j) = data.turbulentOmega(node);
                                    wallDistNow(i,j) = data.wallDistance(node);
%                                     mut(i,j) = data.viscosity_turb(node);
%                                     StR(i,j) = data.strain_rate_mag(node);
                                    if obj.trans == 1
                                        gamma(i,j) = data.intermittency(node);
                                        Reth(i,j) = data.momentum_thickness_re(node);
                                    end
                                end
                            end
                    end
                    obj.ro{nb} = ro;
                    obj.u{nb} = u;
                    obj.v{nb} = v;
                    obj.w{nb} = w;
                    obj.Et{nb} = Et;
                    obj.k{nb} = k;
                    obj.omega{nb} = omega;
                    obj.wallDist{nb} = wallDistNow;
                    if obj.trans == 1
                        obj.gamma{nb} = gamma;
                        obj.Reth{nb} = Reth;
                    end
                end
                obj.mut_store = mut;
                obj.StR_store = StR;
                %clear nodes
                obj.getBCs(blk.inlet_blocks{1});
            end
        end

%         function value = get.p(obj)
%             disp('Calculating p')
%             value = cell(1,obj.NB);
%             for nb = 1:obj.NB
%                 value{nb} = (obj.gas.gam -1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
%             end
%         end

%         function value = get.T(obj)
%             disp('Calculating T')
%             value = cell(1,obj.NB);
%             for nb = 1:obj.NB
%                 pnow = (obj.gas.gam -1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
%                 value{nb} = pnow./(obj.ro{nb}*obj.gas.rgas);
%             end
%         end

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

        end

        function value = get.Pr(obj)
            value = cell(1,obj.NB);
            if ~isempty(obj.Pr_store)
                disp('Using stored Pr')
                value = obj.Pr_store;
            elseif ~isempty(obj.mut_store)
                disp('Calculating Pr: using stored mut')
                for nb = 1:obj.NB
                    value{nb} = obj.mut_store{nb}.*obj.StR{nb}.^2;
    %                 [DUDX,DUDY] = gradHO(obj.blk.x{nb},obj.blk.y{nb},obj.u{nb});
    %                 [DVDX,DVDY] = gradHO(obj.blk.x{nb},obj.blk.y{nb},obj.v{nb});
                    %value{nb} = (obj.mut{nb}.*DUDX.^2 + obj.mut{nb}.*(DUDY.^2+DVDX.^2) + obj.mut{nb}.*DVDY.^2);
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
            value = obj.mut_store;
        end
            
        
    end
end