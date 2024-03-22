classdef DNS_cascade < DNS_case
    %DNS_CASCADE Subclass of DNS_case contaning methods and properties
    ... specific to cascade cases

    properties
    end

    properties (Dependent = true)
        %Re_theta_in
    end

    methods
        function obj = DNS_cascade(casename,run)
            %DNS_CHANNEL Construct an instance of this class
            %   Detailed explanation goes here
            args.topology = 1;
            if nargin > 0
                if nargin < 2
                    run = [];
                end
            else
                casename = [];
                run = [];
            end
            obj@DNS_case(casename,run,args);

            obj.pitch = obj.blk.y{1}(1,end) - obj.blk.y{2}(1,1);
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
            end
            obj.blk.aspect = [(obj.blk.viewarea(2)-obj.blk.viewarea(1)) ...
                (obj.blk.viewarea(4)-obj.blk.viewarea(3)) 1];


            obj.blk.aspect(2) = obj.blk.aspect(2)+obj.pitch;
%             obj.blk.viewarea(4) = obj.blk.viewarea(4)+obj.pitch;
            obj.blk.viewarea(4) = obj.blk.viewarea(4)+0.5*obj.pitch;
%             obj.blk.viewarea(3) = obj.blk.viewarea(3)+0.5*obj.pitch;
            obj.blk.n_pitchwise_repeats = 4;
            obj.blk.pitch = obj.pitch;
        end

        function boundaries = getBoundaries(obj)
            boundaries = {};
            b.label = "Inlet";
            b.blocks = [1 2];
            b.patches = [1 1];
            b.type = 4;
            boundaries{end+1} = b;
            
            b.label = "Upper Periodic";
            b.blocks = [1 5 8];
            b.patches = [4 3 4];
            b.type = 12;
            boundaries{end+1} = b;
            
            b.label = "Lower Periodic";
            b.blocks = [2 6 9];
            b.patches = [3 3 3];
            b.type = 8;
            boundaries{end+1} = b;
            
            b.label = "Outlet";
            b.blocks = [8 9];
            b.patches = [2 2];
            b.type = 5;
            boundaries{end+1} = b;
            
            b.label = "Blade";
            b.blocks = [3 4 5 7];
            b.patches = [4 4 4 4];
            b.type = 3;
            boundaries{end+1} = b;
        end


%         function blkNodes = write_fluent_mesh_2d(obj, path)
%             blkNodes = writeCascadeFluentMesh(path, obj.blk, obj.blk.next_block, obj.blk.next_patch, true);
%         end

        function newCase = instantiate(obj)
            newCase = DNS_channel;
        end

        function pi = P02P01(obj)

            obj.meanFlow.blk.outlet_blocks{1} = [8 9];
            obj.meanFlow.p0out / obj.bcs.Poin;

        end

        function mu_ref = set_Re_c(obj, Re)

            rgas = obj.gas.cp * (1 - 1/obj.gas.gam);
            cut = kCut(obj.blk, obj.gas, obj.bcs);
            M1 = M_p0T0V(obj.bcs.Poin, obj.bcs.Toin, obj.bcs.vin);
            T1 = obj.bcs.Toin * T_T0(M1, obj.gas.gam);
            ro0 = obj.bcs.Poin/(rgas * obj.bcs.Toin);
            ro1 = ro0 * ro_ro0(M1, obj.gas.gam);

            mu1 = ro1*obj.bcs.vin*cut.c / Re;
            mu_ref = sutherland_mu_ref(mu1, T1);


        end




    end
end