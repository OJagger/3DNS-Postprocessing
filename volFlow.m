classdef volFlow < handle
    % VOLFLOW Contains a 3D instntaneous flow from DNS case

        
    properties
        NB;
        gas;
        blk;
        bcs;
        gam;
        cp;
        rgas;
        ro;
        u;
        v;
        w;
        Et;
        mut;
        time;
        tau;
        div;
        vort_x;
        vort_y;
        vort_z;
        rmut;
        nk;
        span;
        z;
        Qstore;              % Q Criterion
        casetype;
        flowpath;
        if_rans;
    end

    properties (Dependent = true, Hidden = true)
        T;              % Temperature
        p;              % p stat
        M;              % Mach No
        s;              % Entropy ( cp*log(T/300) - R*log(p/1e5) )
        vel;            % Velocity
        Msurf;          % Surface Mach No
        MSlice;
        Q;
        
    end

    methods
        function obj = volFlow(casedir, blk, gas, bcs, casetype, if_rans)

            if nargin > 0
                blockdims = blk.blockdims;
                obj.blk = blk;
                obj.gas = gas;
                obj.bcs = bcs;
                obj.nk = blockdims(1,3);
                obj.gam = gas.gam;
                obj.cp = gas.cp;
                obj.rgas = obj.cp*(1-1/obj.gam);
                obj.casetype = casetype;
                obj.if_rans = if_rans;
    
                obj.NB = size(blockdims,1);
                for nb = 1:obj.NB
    
                    fprintf('Reading block %d/%d\n',[nb, obj.NB])

                    ro = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
                    ru = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
                    rv = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
                    rw = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
                    Et = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
                    tau_xx = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
                    tau_yy = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
                    tau_zz = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
                    tau_xy = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
                    tau_xz = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));
                    tau_yz = zeros(blockdims(nb,1),blockdims(nb,2),blockdims(nb,3));

                    switch casetype
                        case 'cpu'
                            if exist(fullfile(casedir, 'inst_flo'), 'dir')
                                flopath = fullfile(casedir, 'inst_flo', ['flo2_' num2str(nb)]);
                                flofile = fopen(flopath,'r');
                                nodfile = fopen(fullfile(casedir, 'inst_flo', ['nod2_' num2str(nb)]),'r');
                            else
                                flopath = fullfile(casedir,  ['flo2_' num2str(nb)]);
                                flofile = fopen(flopath,'r');
                                nodfile = fopen(fullfile(casedir, ['nod2_' num2str(nb)]),'r');
                            end
                            %viscpath = fullfile(casedir,  ['visc_' num2str(nb)]);
                            %viscfile = fopen(viscpath,'r');
        
                            A = fread(flofile,inf,'float64');
                            A = reshape(A,5,length(A)/5);
                            
                            B = fread(nodfile,inf,'uint32');
                            B = reshape(B,3,length(B)/3);
        
        %                     C = fread(viscfile,inf,'float64');
        %                     C = reshape(C,13,length(C)/13 );
        
                    
                            fclose(flofile);
                            fclose(nodfile);
            
                            for n=1:size(A,2)
                                i = B(1,n);
                                j = B(2,n);
                                k = B(3,n);
                                ro(i,j,k) = A(1,n);
                                ru(i,j,k) = A(2,n);
                                rv(i,j,k) = A(3,n);
                                rw(i,j,k) = A(4,n);
                                Et(i,j,k) = A(5,n);
                                mut = [];
                            end
                        case 'gpu'
                            ni = blk.blockdims(nb,1);
                            nj = blk.blockdims(nb,2);
                            nk = blk.blockdims(nb,3);
                            fid = fopen(fullfile(casedir, ['flow_' num2str(nb)]));
                            A = fread(fid, ni*nj*nk*5, 'float64');
                            A = reshape(A, 5, length(A)/5)';
                            fclose(fid);

                            ro = reshape(A(:,1),ni,nj,nk);
                            ru = reshape(A(:,2),ni,nj,nk);
                            rv = reshape(A(:,3),ni,nj,nk);
                            rw = reshape(A(:,4),ni,nj,nk);
                            Et = reshape(A(:,5),ni,nj,nk);

                            if if_rans
                                fid = fopen(fullfile(casedir, ['rans_' num2str(nb)]));
                                A = fread(fid, ni*nj*nk, 'float64');
                                mut = reshape(A(:),ni,nj,nk);
                                fclose(fid);
                            else
                                mut = [];
                            end
                    end
    
                    obj.ro{nb} = ro;
                    obj.u{nb} = ru./ro;
                    obj.v{nb} = rv./ro;
                    obj.w{nb} = rw./ro;
                    obj.Et{nb} = Et;
                    obj.mut{nb} = mut;
                end
            else
                obj.ro = {};
                obj.u = {};
                obj.v = {};
                obj.w = {};
                obj.Et = {};
                obj.mut = {};
            end
        end

        function newFlow = copySkeleton(obj)
            newFlow = volFlow;
            newFlow.NB = obj.NB;
            newFlow.gas = obj.gas;
            newFlow.blk = obj.blk;
            newFlow.gam = obj.gam;
            newFlow.cp = obj.cp;
            newFlow.rgas = obj.rgas;
            newFlow.time = obj.time;
        end

        function value = get.p(obj)
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                value{nb} = (obj.gam -1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
            end
        end

        function value = get.T(obj)
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                value{nb} = obj.p{nb}./(obj.ro{nb}*obj.rgas);
            end
        end

        function value = get.vel(obj)
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                value{nb} = sqrt(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2);
            end
        end

        function value = get.M(obj)
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                value{nb} = obj.vel{nb}./sqrt(obj.gam*obj.rgas*obj.T{nb});
            end
        end

        function value = get.s(obj)
            pool = gcp('nocreate');
            value = cell(1,obj.NB);
            if isempty(pool)
                for nb =1:obj.NB
                    pnow = (obj.gam -1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                    Tnow = pnow./(obj.ro{nb}*obj.rgas);
                    value{nb} = obj.cp*log(Tnow/300) - obj.rgas*log(pnow/1e5);
                end
            else
                parfor nb =1:obj.NB
                    pnow = (obj.gam -1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb})
                    Tnow = pnow./(obj.ro{nb}*obj.rgas);
                    value{nb} = obj.cp*log(Tnow/300) - obj.rgas*log(pnow/1e5);
                end
            end
        end

        function value = get.Msurf(obj)
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                value{nb} = obj.cp*log(obj.T{nb}/300) - obj.rgas*log(obj.p{nb}/1e5);
            end
        end

        function value = get.Q(obj)
            value = obj.get_Q;
        end

        function obj = set.Q(obj, value)
            obj.set_Q(value);
        end

        function set_Q(obj, value)
        end

        function value = get_Q(obj)
            if isempty(obj.Qstore)
                value = cell(1,obj.NB);
                for nb = obj.blk.oblocks
                    fprintf('Calculating Q criterion in block %d\n',nb)
                    value{nb} = Q_criterion(obj.blk.x{nb},obj.blk.y{nb},obj.blk.z,obj.u{nb},obj.v{nb},obj.w{nb});
                end
                obj.Qstore = value;
            else
                value = obj.Qstore;
            end
        end

        function s = plot_Q_criterion(obj, thresh)
%             X = repmat(obj.blk.x{1},[1 1 obj.nk]);
%             Y = repmat(obj.blk.y{1},[1 1 obj.nk]);
%             Z = repmat(obj.blk.z,[size(obj.blk.x{1}) 1]);
            Qnow = obj.Q;
            for nb = obj.blk.oblocks
                %X = repmat(obj.blk.x{nb},1,1,length(obj.blk.z));
                X = linspace(0,1,size(obj.blk.x{nb},1));
                %Y = repmat(obj.blk.y{nb},1,1,length(obj.blk.z));
                Y = linspace(0,1,size(obj.blk.y{nb},2));
                %znow = reshape(obj.blk.z,1,1,[]);
                %Z = repmat(znow, size(X,1), size(X,2));
                Z = linspace(0,1,length(obj.blk.z));
                %[X, Y, Z] = meshgrid(X,Y,Z);
                Qblock = permute(Qnow{nb},[2 1 3]);
                s = isosurface(X, Y, Z, Qblock, thresh);
                %p = patch(s);
            end
            %axis equal
        end

        function value = k_ave(obj,prop)
            pool = gcp('nocreate');
            if isempty(pool)
                for nb = 1:obj.NB
                    switch prop
                        case 's'
                            pnow = (obj.gam -1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                            Tnow = pnow./(obj.ro{nb}*obj.rgas);
                            propnow = obj.cp*log(Tnow/300) - obj.rgas*log(pnow/1e5);
                    end
                    value{nb} = mean(propnow,3);
                end
            else
                parfor nb = 1:obj.NB
                    switch prop
                        case 's'
                            pnow = (obj.gam -1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                            Tnow = pnow./(obj.ro{nb}*obj.rgas);
                            propnow = obj.cp*log(Tnow/300) - obj.rgas*log(pnow/1e5);
                    end
                    value{nb} = mean(propnow,3);
                end
            end
        end

            

        function getSize(obj)
            props = properties(obj); 
            totSize = 0; 
           
            for ii=1:length(props) 
                currentProperty = getfield(obj, char(props(ii))); 
                temp = whos('currentProperty'); 
                totSize = totSize + temp.bytes; 
            end
          
            fprintf(1, '%d MB\n', totSize/1e6);
        end

        function value = slice(obj,k)
            if nargin == 1
                k = floor(obj.nk/2);
            end

            if obj.if_rans
                value = RANSSlice(obj.blk,obj.gas,obj.bcs);
            else
                value = kSlice(obj.blk, obj.gas, obj.bcs); %(obj.blk.blockdims, obj.gas);
            end

            for nb=1:obj.NB
                nb;
                value.ro{nb} = obj.ro{nb}(:,:,k);
                value.u{nb} = obj.u{nb}(:,:,k);
                value.v{nb} = obj.v{nb}(:,:,k);
                value.w{nb} = obj.w{nb}(:,:,k);
                value.Et{nb} = obj.Et{nb}(:,:,k);
                if obj.if_rans
                    value.mut_store{nb} = obj.mut{nb};
                end
            end
        end

        function write_tecplot_files(obj, path)
            for nb = 1:obj.NB
                    tdata=[];
                    tdata.Nvar=9;
                    tdata.varnames={'x','y','z','p','T','ro','u','v','w'};
                    tdata.cubes(1).zonename=['block ',num2str(nb)];
                    zz(1,1,1:obj.blk.nk) = obj.blk.z; 
                    z = repmat(zz,[obj.blk.blockdims(nb,1) obj.blk.blockdims(nb,2) 1]);
                    x = repmat(obj.blk.x{nb},[1 1 obj.blk.nk]);
                    y = repmat(obj.blk.y{nb},[1 1 obj.blk.nk]);
                if ismember(nb, obj.blk.oblocks_flip)
                    tdata.cubes(1).x=x(end:-1:1,:,:);
                    tdata.cubes(1).y=y(end:-1:1,:,:);
                    tdata.cubes(1).z=z(end:-1:1,:,:);
                    tdata.cubes(1).v(1,:,:,:)=obj.p{nb}(end:-1:1,:,:);
                    tdata.cubes(1).v(2,:,:,:)=obj.T{nb}(end:-1:1,:,:);
                    tdata.cubes(1).v(3,:,:,:)=obj.ro{nb}(end:-1:1,:,:);
                    tdata.cubes(1).v(4,:,:,:)=obj.u{nb}(end:-1:1,:,:);
                    tdata.cubes(1).v(5,:,:,:)=obj.v{nb}(end:-1:1,:,:);
                    tdata.cubes(1).v(6,:,:,:)=obj.w{nb}(end:-1:1,:,:);
                else
                    tdata.cubes(1).x=x;
                    tdata.cubes(1).y=y;
                    tdata.cubes(1).z=z;
                    tdata.cubes(1).v(1,:,:,:)=obj.p{nb};
                    tdata.cubes(1).v(2,:,:,:)=obj.T{nb};
                    tdata.cubes(1).v(3,:,:,:)=obj.ro{nb};
                    tdata.cubes(1).v(4,:,:,:)=obj.u{nb};
                    tdata.cubes(1).v(5,:,:,:)=obj.v{nb};
                    tdata.cubes(1).v(6,:,:,:)=obj.w{nb};
                end
                tdata.vformat(1:9) = 2; 
                tdata.cubes.solutiontime=0;
                fname = fullfile(path,['tec_flow_',num2str(nb),'.plt']);
                mat2tecplot(tdata,fname)
            end
        end

        function interpOntoNewGrid(obj, newcase)

            newFlow = volFlow();
            ronow = {};


            for ib = 1:obj.NB
                fprintf('interpolting block: %d\n', ib)

                [fic, fjc] = obj.get_spacing(obj.blk.x{ib},obj.blk.y{ib});
                [fi, fj] = obj.get_spacing(newcase.blk.x{ib},newcase.blk.y{ib});
                fkc = linspace(0,1,obj.blk.nk);
                fk = linspace(0,1,newcase.blk.nk);

                [Jc,Ic,Kc] = meshgrid(fjc,fic,fkc);
                [J,I,K] = meshgrid(fj,fi,fk);

                ronow{ib} = interp3(Jc,Ic,Kc,obj.ro{ib},J,I,K);
                unow{ib} = interp3(Jc,Ic,Kc,obj.u{ib},J,I,K);
                vnow{ib} = interp3(Jc,Ic,Kc,obj.v{ib},J,I,K);
                wnow{ib} = interp3(Jc,Ic,Kc,obj.w{ib},J,I,K);
                Etnow{ib} = interp3(Jc,Ic,Kc,obj.Et{ib},J,I,K);
                
            end

            newFlow = obj;
            newFlow.ro = ronow;
            newFlow.u = unow;
            newFlow.v = vnow;
            newFlow.w = wnow;
            newFlow.Et = Etnow;

            newcase.instFlow = newFlow;
%             newcase.ins
        end


        function writeFlowBlock(obj, path, nb, type)
            
            ronow = obj.ro{nb};
            runow = ronow.*obj.u{nb};
            rvnow = ronow.*obj.v{nb};
            rwnow = ronow.*obj.w{nb};
            Etnow = obj.Et{nb};

            switch obj.casetype
                case 'cpu'
                    fpos = 0;
                    for i=1:size(obj.ro{nb},1)
                        for j=1:size(obj.ro{nb},2)
                            for k=1:size(obj.ro{nb},3)
                                fpos = fpos+1;
                                B(1, fpos) = i;
                                B(2, fpos) = j;
                                B(3,fpos) = k;
        
                                A(1,fpos) = ronow(i,j,k);
                                A(2,fpos) = runow(i,j,k);
                                A(3,fpos) = rvnow(i,j,k);
                                A(4,fpos) = rwnow(i,j,k);
                                A(5,fpos) = Etnow(i,j,k);
                            end
                        end
                    end
        
                    flopath = fullfile(path,  ['flo2_' num2str(nb)]);
                    flofile = fopen(flopath,'w');
                    nodfile = fopen(fullfile(path, ['nod2_' num2str(nb)]),'w');
                    %viscpath = fullfile(casedir,  ['visc_' num2str(nb)]);
                    %viscfile = fopen(viscpath,'r');
        
                    fwrite(flofile,A,'float64');
                    fwrite(nodfile,B,'uint32');
        
            
                    fclose(flofile);
                    fclose(nodfile);
                case 'gpu'

                            A(:,:,:,1) = ronow;
                            A(:,:,:,2) = runow;
                            A(:,:,:,3) = rvnow;
                            A(:,:,:,4) = rwnow;
                            A(:,:,:,5) = Etnow;
                            A = permute(A,[4 1 2 3]);
                            

                            fid = fopen(fullfile(path, ['flow_' num2str(nb)]),'w');
                            fwrite(fid,A,'float64');
                            fclose(fid)
                    end
        end

        function writeFlow(obj, path)
            if nargin < 2
                path = obj.flowpath;
            end
            for nb = 1:obj.NB
                fprintf('Writing flow in block %d/%d\n',[nb,obj.NB])
                obj.writeFlowBlock(path, nb);
            end
        end 

        function newFlow = sub_sample_block(obj, block, varargin)
            p = inputParser;
            addRequired(p,'block');
            addOptional(p,'factor',1);
            addOptional(p,'dims',[1 2 3]);

            newFlow = obj.copySkeleton;
            newFlow.NB = 1;
            parse(p, block, varargin{:});
            ib = p.Results.block;
            factor = p.Results.factor;
            
            if factor ~= 1
                is = 1:factor:obj.blk.blockdims(ib,1);
                js = 1:factor:obj.blk.blockdims(ib,2);
                ks = 1:factor:obj.blk.blockdims(ib,3);
            else
                is = 1:dims(1):obj.blk.blockdims(ib,1);
                js = 1:dims(2):obj.blk.blockdims(ib,2);
                ks = 1:dims(3):obj.blk.blockdims(ib,3);
            end

            blk = obj.blk;
            blk.blockdims = [length(is) length(js) length(ks)];
            blk.x{1:obj.NB} = []; blk.x{ib} = obj.blk.x{ib}(is,js);
            blk.y{1:obj.NB} = []; blk.y{ib} = obj.blk.y{ib}(is,js);
            blk.nk = length(ks);
            blk.z = blk.z(ks);

            newFlow.blk = blk;
            newFlow.nk = blk.nk;
            newFlow.ro{ib} = obj.ro{ib}(is,js,ks);
            newFlow.u{ib} = obj.u{ib}(is,js,ks);
            newFlow.v{ib} = obj.v{ib}(is,js,ks);
            newFlow.w{ib} = obj.w{ib}(is,js,ks);
            newFlow.Et{ib} = obj.Et{ib}(is,js,ks);

        end
        
    end
    methods (Static)
        function [fi, fj] = get_spacing(x,y)
            iq = ceil(size(x,1)/2);
            xj = x(iq,:);
            yj = y(iq,:);

            jq = ceil(size(x,2)/2);
            xi = x(:,jq);
            yi = y(:,jq);

            si = path_length(xi,yi);
            sj = path_length(xj,yj);

            fi = si/si(end);
            fj = sj/sj(end);
        end
    end

end