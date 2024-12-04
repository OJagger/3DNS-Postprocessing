classdef jSlice < handle
    % JSLICE Contains a 2D slice of the flow next to the blade surface
    %   Detailed explanation goes here

    properties
        NB;
        blk;
        gas;
        bcs;
        ro;
        u;
        v;
        w;
        Et;
        time;
        nSlice;
        ssurf;
        Y0;
        X;
        Z;
        pdyn;              % Freestream dynamic pressure
        casepath;
    end

    properties (Dependent = true)
        T;
        p;
        M;
        s;
        vel;
        mu;
        tau_w;
        cf;
    end

    methods
        function obj = jSlice(blk, gas, bcs, casedir, nSlice, time, casetype, ishere)
            
            if nargin < 8
                ishere = false;
            end

            if nargin > 0

                obj.blk = blk;
                obj.gas = gas;
                obj.bcs = bcs;

            end

            if nargin > 3
                obj.nSlice = nSlice;
                obj.time = time;
                obj.ro = [];
                obj.u = [];
                obj.v = [];
                obj.w = [];
                obj.Et = [];
                obj.NB = size(blk.blockdims,1);
                x = [];
                y = [];
                y0 = [];
                if ~ishere
                    casedir = [casedir '/j_cuts'];
                end
                for iblk=1:length(blk.oblocks)
                    nb = blk.oblocks(iblk);
                    flopath = fullfile(casedir, ['jcu2_' num2str(nb) '_' num2str(nSlice)]);
                    nodpath = fullfile(casedir, ['jnd2_' num2str(nb) '_' num2str(nSlice)]);
                    flofile = fopen(flopath,'r');
                    nodfile = fopen(nodpath,'r');
                    A = fread(flofile,inf,'float64');
                    A = reshape(A,5,length(A)/5);
                    
                    B = fread(nodfile,inf,'uint32');
                    B = reshape(B,3,length(B)/3);
            
                    fclose(flofile);
                    fclose(nodfile);
    
                    ro = zeros(blk.blockdims(nb,1),blk.blockdims(nb,3));
                    ru = zeros(blk.blockdims(nb,1),blk.blockdims(nb,3));
                    rv = zeros(blk.blockdims(nb,1),blk.blockdims(nb,3));
                    rw = zeros(blk.blockdims(nb,1),blk.blockdims(nb,3));
                    Et = zeros(blk.blockdims(nb,1),blk.blockdims(nb,3));
                    
    
                    for n=1:size(A,2)
                        i = B(1,n);
                        j = B(2,n);
                        k = B(3,n);
                        ro(i,k) = A(1,n);
                        ru(i,k) = A(2,n);
                        rv(i,k) = A(3,n);
                        rw(i,k) = A(4,n);
                        Et(i,k) = A(5,n);
                    end

                    xtmp = blk.x{nb}(:,end-1);
                    ytmp = blk.y{nb}(:,end-1);
                    dxtmp = blk.x{nb}(:,end-1) - blk.x{nb}(:,end);
                    dytmp = blk.y{nb}(:,end-1) - blk.y{nb}(:,end);
                    y0tmp = sqrt(dxtmp.^2+dytmp.^2);
                    
                    if blk.oblocks_flip(iblk) == 1
                        ro = flip(ro);
                        ru = flip(ru);
                        rv = flip(rv);
                        rw = flip(rw);
                        Et = flip(Et);
                        xtmp = flip(xtmp);
                        ytmp = flip(ytmp);
                        y0tmp = flip(y0tmp);
                    end
                    
                    obj.ro = [obj.ro; ro(1:end-1,:)];
                    obj.u = [obj.u; ru(1:end-1,:)./ro(1:end-1,:)];
                    obj.v = [obj.v; rv(1:end-1,:)./ro(1:end-1,:)];
                    obj.w = [obj.w; rw(1:end-1,:)./ro(1:end-1,:)];
                    obj.Et = [obj.Et; Et(1:end-1,:)];
                    x = [x; xtmp(1:end-1)];
                    y = [y; ytmp(1:end-1)];
                    y0 = [y0; y0tmp];

                end
                x = x';
                y = y';
                y0 = y0';
                [~, iLE] = min(x);
                [~, iTE] = max(x);
                obj.ro = obj.ro(iLE:iTE,:);
                obj.u = obj.u(iLE:iTE,:);
                obj.v = obj.v(iLE:iTE,:);
                obj.w = obj.w(iLE:iTE,:);
                obj.Et = obj.Et(iLE:iTE,:);

                x = x(iLE:iTE);
                y = y(iLE:iTE);
                y0 = y0(iLE:iTE);
                obj.ssurf(1) = 0;
                for i=2:length(x)
                    dx = x(i)-x(i-1);
                    dy = y(i)-y(i-1);
                    ds = sqrt(dx^2 + dy^2);
                    obj.ssurf(i) = obj.ssurf(i-1)+ds;
                end
                z = linspace(0,blk.span,blk.nk);
                [obj.Z, obj.X] = meshgrid(z,obj.ssurf);
                obj.Y0 = repmat(y0',1,blk.nk);
            end
        end

        function value = get.p(obj)
            value = (obj.gas.gam -1)*(obj.Et - 0.5*(obj.u.^2 + obj.v.^2 + obj.w.^2).*obj.ro);
        end

        function value = get.T(obj)
            value = obj.p./(obj.ro*obj.gas.rgas);
        end

        function value = get.vel(obj)
            value = sqrt(obj.u.^2 + obj.v.^2 + obj.w.^2);
        end

        function value = get.M(obj)
            value  = obj.vel./sqrt(obj.gas.gam*obj.gas.rgas*obj.T);
        end

        function value = get.s(obj)
            value = obj.gas.cp*log(obj.T/300) - obj.gas.rgas*log(obj.p/1e5);
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

        function jPlot(obj,prop,ax,lims,label)
            if nargin < 3 || isempty(ax)
                ax = gca;
            end
            q = obj.(prop);
            dx = abs(obj.X(end,1)-obj.X(1,1));
            dy = abs(obj.Z(1,end)-obj.Z(1,1));
            
            pcolor(ax, obj.X, obj.Z, q);
            xlabel('Surface distance')
            shading('interp')
            axis equal
            
            cb = colorbar('southoutside');
            if nargin > 3 && ~isempty(lims)
                caxis(lims);
            end
            if nargin > 4 && ~isempty(label)
                cb.Label.String = label;
                cb.FontSize = 9;
            end
            pbaspect(ax, [dx dy dx]);
            set(ax, 'FontSize', 9)
        end

        function value = get.mu(obj)
            disp('Calcualting mu')
            pnow = (obj.gas.gam - 1)*(obj.Et - 0.5*(obj.u.^2 + obj.v.^2 + obj.w.^2).*obj.ro);
            Tnow = pnow./(obj.ro*obj.gas.rgas);
            value = obj.gas.mu_ref*(Tnow/obj.gas.mu_tref).^(3/2) .* (obj.gas.mu_cref + obj.gas.mu_tref)./(obj.gas.mu_cref + Tnow);
        end

        function value = get.tau_w(obj)
            value = obj.mu.*obj.vel./obj.Y0;
        end

        function value = get.cf(obj)
            q = obj.pdyn;
            if ~isempty(q)
                value = obj.tau_w./repmat(obj.pdyn',[1 size(obj.tau_w,2)]);
            else
                fprintf('Set jSlice dynamic pressure before calculating cf\n')
            end
        end
    end
end