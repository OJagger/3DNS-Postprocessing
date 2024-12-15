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

        function Re = Re_c(obj)

            cut = kCut(obj.blk, obj.gas, obj.bcs);
            M1 = M_p0T0V(obj.bcs.Poin, obj.bcs.Toin, obj.bcs.vin);
            T1 = obj.bcs.Toin * T_T0(M1, obj.gas.gam);
            ro0 = obj.bcs.Poin/(obj.gas.rgas * obj.bcs.Toin);
            ro1 = ro0 * ro_ro0(M1, obj.gas.gam);

            mu1 = sutherland_mu(T1, obj.gas.mu_ref, obj.gas.mu_cref, obj.gas.mu_tref);
            Re = ro1*obj.bcs.vin*cut.c/mu1;

        end

        function [y, prof] = wake_profile(obj, slice, xs, prop)

            ymin = 1e4;
            ymax = -1e4;

            for ib = 1:obj.NB
                ymin = min(ymin, min(obj.blk.y{ib}, [], 'all'));
                ymax = max(ymax, max(obj.blk.y{ib}, [], 'all'));
            end

            n = 401;
            x(1:n) = xs;
            y=linspace(ymin, ymax, n);

            ys = slice.unstructured_sample(x, y, 'y');
            y = linspace(min(ys), max(ys), n);

            prof = slice.unstructured_sample(x, y, prop);
        end

        function sol = write_MISES_input(obj)

            f = fopen(fullfile(obj.casdir, 'stream.mises'));


            f = fopen(fullfile(obj.casdir, 'ises.mises'));

        end

        function [bl1, bl2] = read_MISES_bl(obj,ext,path)
            
            if nargin< 3
                fpath = fullfile(obj.casepath, 'MISES', ext);
            else
                fpath = path;
            end

            bl1 = read_mrchbl_output(fullfile(fpath, 'bl1.mises'));
            bl1.s = bl1.x;
            Idat = mis_read_idat('mises',fpath);
            Ises = mis_read_ises('mises',fpath);

            xstag = interp1(Idat.sb{1}, Idat.xb{1}, Idat.sble);
            ystag = interp1(Idat.sb{1}, Idat.yb{1}, Idat.sble);
            [~, istag] = max(Idat.sb{1}.*(Idat.sb{1}<Idat.sble));
            xss = [xstag Idat.xb{1}(istag:-1:1)];
            yss = [ystag Idat.yb{1}(istag:-1:1)];
            sb = Idat.sble - Idat.sb{1};
            sb = [0 sb(istag:-1:1)];
            scale = obj.pitch/Idat.pitch;
            bl1.pitch = Idat.pitch;

            cut = kCut(obj.blk, obj.gas, obj.bcs);

            xb = interp1(sb,xss,bl1.x);
            yb = interp1(sb,yss,bl1.x);

            bl1.s = bl1.s*scale;
            [xle, ile] = min(xb);
            yle = yb(ile);
            bl1.x = scale*(xb-xle);
            bl1.y = scale*(yb-yle);
            bl1.delStar = scale*bl1.delStar;
            bl1.theta = scale*bl1.theta;
            bl1.cax = 1/scale;
            bl1.c = cut.c/scale;

            bl1.Res = Ises.reyn*bl1.s;
            bl1.Re_c = Ises.reyn*bl1.c;

            if isfield(bl1,'xShock')
                bl1.sShock = bl1.xShock;
                bl1.xShock = interp1(sb,xss,bl1.sShock);
            end

            bl2 = read_mrchbl_output(fullfile(fpath, 'bl2.mises'));
            bl2.s = bl2.x;
            Idat = mis_read_idat('mises',fpath);
            Ises = mis_read_ises('mises',fpath);

            xps = [xstag Idat.xb{1}(istag+1:end)];
            yps = [ystag Idat.yb{1}(istag+1:end)];
            sb = Idat.sb{1} - Idat.sble;
            sb = [0 sb(istag+1:end)];
            scale = obj.pitch/Idat.pitch;
            bl2.pitch = Idat.pitch;


            xb = interp1(sb,xps,bl2.x);
            yb = interp1(sb,yps,bl2.x);

            bl2.s = bl2.s*scale;
            bl2.x = scale*(xb-xle);
            bl2.y = scale*(yb-yle);
            bl2.delStar = scale*bl2.delStar;
            bl2.theta = scale*bl2.theta;
            bl2.cax = 1/scale;
            bl2.c = cut.c/scale;

            bl2.Res = Ises.reyn*bl2.s;
            bl2.Re_c = Ises.reyn*bl2.c;

        end
            
    end
end