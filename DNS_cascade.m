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
            
        end


        function blkNodes = write_fluent_mesh_2d(obj, path)
            blkNodes = writeCascadeFluentMesh(path, obj.blk, obj.blk.next_block, obj.blk.next_patch, true);
        end

        function newCase = instantiate(obj)
            newCase = DNS_channel;
        end



    end
end