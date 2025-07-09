classdef kCut < flowSlice
    %KCUT Superclass for instantaneous and mean flow on a k boundary

    properties
        xSurf;
        xBL;
        yBL;
        xO;
        yO;
        iO;
        jO;
        xLE;
        yLE;
        xTE;
        yTE;
        blkO;
        oblocks;
        oblocks_flip;
        iLE;
        iTE;
        niBL;
        c;
        n;
        ssurf;          % Surface distance fron LE
        owrapblock = 0; % wall topology - 0 for channel, nb of last block in o
                        % grid for closed loop (eg blade profile)
        iSS = false;
    end

    properties (Dependent = true, Hidden = true)
        vortZ;          % Z vorticity
        aspect_ratio;
        DUDX;
        DUDY;
        DVDX;
        DVDY;
    end

    methods
        function obj = kCut(blk, gas, bcs, iPS)
            obj@flowSlice(blk,gas,bcs)
            disp('Constructing kCut')
            if nargin > 0

                if nargin > 3
                    if islogical(iPS)
                        obj.iSS = ~iPS;
                    elseif any(strcmp(iPS,{'ps', 'PS'}))
                        obj.iSS = false;
                    else
                        obj.iSS = true;
                    end
                else
                    obj.iSS = true;
                end
                
                if size(blk.blockdims,1) > 0
                    obj.oblocks = blk.oblocks;
                    obj.oblocks_flip = blk.oblocks_flip;
                    obj.gas = gas;
                    obj.blk = blk;
                    %obj.gas.mu_ref = 0.0008748708280693193;
                    obj.gas.rgas = obj.gas.cp*(1-1/obj.gas.gam);
                    obj.NB = size(blk.blockdims,1);
                    xo = [];
                    yo = [];
                    io = [];
                    jo = [];
                    blko = [];
                    for i=1:length(obj.oblocks)
                        nb = obj.oblocks(i);
                        ni = size(blk.x{nb},1);
                        nj = size(blk.x{nb},2);
                        xtmp = blk.x{nb}(:,:);
                        ytmp = blk.y{nb}(:,:);
                        itmp = 1:ni;
                        itmp = repmat(itmp',[1 nj]);
                        jtmp = 1:nj;
                        if blk.next_patch{nb}.jp == 3
                            jtmp = flip(jtmp);
                        end
                        jtmp = repmat(jtmp, [ni 1]);
                        blktmp = nb*ones(ni,nj);
                        if blk.next_patch{nb}.jp == 3
                            xtmp = flip(xtmp,2);
                            ytmp = flip(ytmp,2);
                        end
                        if obj.oblocks_flip(i) == 1
                            xtmp = flip(xtmp);
                            ytmp = flip(ytmp);
                            itmp = flip(itmp);
                        end
                        if size(xo,1) == 0
                            xo = xtmp; yo = ytmp;
                            io = itmp; jo = jtmp;
                            blko = blktmp;
                        else
                            xo = [xo; xtmp(2:end,:)];
                            yo = [yo; ytmp(2:end,:)];
                            io = [io; itmp(2:end,:)];
                            jo = [jo; jtmp(2:end,:)];
                            blko = [blko; blktmp(2:end,:)];
                        end
                        
                        if xo(1,1) == xo(end,1)
                            xo = xo(1:end-1,:);
                            yo = yo(1:end-1,:);
                            io = io(1:end-1,:);
                            jo = jo(1:end-1,:);
                            blko = blko(1:end-1,:);
                            obj.owrapblock = nb;
                        end

                    end
                    xsurf = xo(:,1);
                    [~, obj.iLE] = min(xsurf);
                    [~, obj.iTE] = max(xsurf);
                    if obj.iSS
                        obj.xSurf = xsurf(obj.iLE:obj.iTE);
                        obj.xO = xo(obj.iLE:obj.iTE,:);
                        obj.yO = yo(obj.iLE:obj.iTE,:);
                        obj.iO = io(obj.iLE:obj.iTE,:);
                        obj.jO = jo(obj.iLE:obj.iTE,:);
                        obj.blkO = blko(obj.iLE:obj.iTE,:);
                   else
                        obj.xSurf = xsurf([obj.iLE:-1:1 end:-1:obj.iTE]);
                        obj.xO = xo([obj.iLE:-1:1 end:-1:obj.iTE],:);
                        obj.yO = yo([obj.iLE:-1:1 end:-1:obj.iTE],:);
                        obj.iO = io([obj.iLE:-1:1 end:-1:obj.iTE],:);
                        obj.jO = jo([obj.iLE:-1:1 end:-1:obj.iTE],:);
                        obj.blkO = blko([obj.iLE:-1:1 end:-1:obj.iTE],:);
                    end
                    obj.c = sqrt((obj.xO(end,1) - obj.xO(1,1))^2 + (obj.yO(end,1) - obj.yO(1,1))^2);
                    %obj.xSurf = xsurf([obj.iLE:-1:1 end:-1:obj.iTE]);
                    %size(obj.xSurf)
                    R = [0 -1; 1 0];
                    obj.xLE = obj.xO(1,1);
                    obj.yLE = obj.yO(1,1);
                    obj.xTE = obj.xO(end,1);
                    obj.yTE = obj.yO(end,1);

                    if obj.iSS
                        isurf = obj.iLE:obj.iTE;
                    else
                        isurf = [obj.iLE:-1:1 size(xo,1):-1:obj.iTE];
                    end
                    obj.yBL = zeros(size(obj.xO));
                    obj.ssurf = zeros(1,size(obj.yBL,1));
                    ii = 0;
                    for i = isurf
                        ii = ii +1;
                        if i ==1
                            s1 = [(xo(i+1,1)-xo(i,1)); (yo(i+1,1)-yo(i,1))];
                            n1 = R*s1/norm(s1);
                            n2 = n1;
                        elseif i == size(xo,1)
                            s2 = [(xo(i,1)-xo(i-1,1)); (yo(i,1)-yo(i-1,1))];
                            n2 = R*s2/norm(s2);
                            n1 = n2;
                        else
                            s1 = [(xo(i+1,1)-xo(i,1)); (yo(i+1,1)-yo(i,1))];
                            s2 = [(xo(i,1)-xo(i-1,1)); (yo(i,1)-yo(i-1,1))];
                            n1 = R*s1/norm(s1);
                            n2 = R*s2/norm(s2);
                        end
                        nnow = (0.5*(n1+n2)/norm(0.5*(n1+n2)));
                        obj.n(:,ii) = nnow;
                        for j=1:size(xo,2)
                            dx = xo(i,j) - xo(i,1);
                            dy = yo(i,j) - yo(i,1);
                            obj.yBL(ii,j) = dot(nnow, [dx;dy]);
                            obj.xBL(ii,j) = obj.xSurf(ii);
                        end
                        if ii>1
                            dx = obj.xO(ii,1) - obj.xO(ii-1,1);
                            dy = obj.yO(ii,1) - obj.yO(ii-1,1);
                            ds = sqrt(dx^2 + dy^2);
                            obj.ssurf(ii) = obj.ssurf(ii-1) + ds;
                        end
                    end
                    obj.niBL = length(obj.xSurf);
                else
                    obj.gas = gas;
                    obj.gas.rgas = obj.gas.cp*(1-1/obj.gas.gam);
                    obj.NB = size(blk.blockdims,1);
                    obj.blk = blk;
                end

            end
        end             % End of constructor

        function value = get.vortZ(obj)
            disp('Calculating z componant of vorticity')
            value = cell(1,obj.NB);
            for nb = 1:obj.NB
                [~, dudy] = gradHO(obj.blk.x{nb},obj.blk.y{nb},obj.u{nb});
                [dvdx, ~] = gradHO(obj.blk.x{nb},obj.blk.y{nb},obj.v{nb});
                value{nb} = dvdx - dudy;
            end
        end

        function value = get.DUDX(obj)
            for ib=1:obj.NB
                [value{ib},~] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.u{ib});
            end
        end

        function value = get.DUDY(obj)
            for ib=1:obj.NB
                [~,value{ib}] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.u{ib});
            end
        end

        function value = get.DVDX(obj)
            for ib=1:obj.NB
                [value{ib},~] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.v{ib});
            end
        end

        function value = get.DVDY(obj)
            for ib=1:obj.NB
                [~,value{ib}] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.v{ib});
            end
        end


        function [profile, i] = BLprof(obj, x, prop)
            %BLPROF get a profile of prop across the BL at specified x
            i = obj.x2ind(x);
            if any(strcmp(["dsdy" "U" "yBL" "yplus"],prop))
                profile = obj.(prop)(i,:);
                size(profile)
            else
                propfield = oGridProp(obj,prop);
                profile = propfield(i,:);
            end
        end

        function value = surface_normal(obj, x)
            i = obj.x2ind(x);
            
            value = obj.n(:,i);
        end

        function blfield = oGridProp(obj, prop)
            %OGRIDPROP Construct array of property in o grid for one surface
            %of blade

            propnow = obj.(prop);
            if ~iscell(propnow)
            
            %if any(strcmp(["U","dsdy","yBL","yplus"], prop))
                blfield = propnow;
            else
                blfield = [];
                propnow = obj.(prop);
                for i=1:length(obj.oblocks)
                    nb = obj.oblocks(i);
                    clear temp
                    temp = propnow{nb};
                    if obj.blk.next_patch{nb}.jp == 3
                        temp = flip(temp,2);
                    end
                    %size(temp)
                    if obj.oblocks_flip(i) == 1
                        temp = flip(temp);
                    end
                    if size(blfield,1) == 0
                        blfield = temp;
                    elseif nb == obj.owrapblock
                        blfield = [blfield; temp(2:end-1,:,:,:)];
                    else
                        blfield = [blfield; temp(2:end,:,:,:)];
                    end
                end
                if obj.iSS
                    blfield = blfield(obj.iLE:obj.iTE,:,:,:);
                else
                    blfield = blfield([obj.iLE:-1:1 end:-1:obj.iTE],:,:,:);
                end
                %size(blfield,1)
                %blfield = blfield([obj.iLE:-1:1 end:-1:obj.iTE],:);
                %size(blfield,1)
            end
        end

        function value = hGridProp(obj, prop)
            value = cell(1,obj.NB);
            v = obj.(prop);
            if size(v,2) == 1 || size(v,1) == 1
                for i=1:length(v)
                    value{obj.blkO(i,1)}(obj.iO(i,1)) = v(i);
                end
            else
                for i=1:size(v,1)
                    for j=1:size(v,2)
                        value{obj.blkO(i,1)}(obj.iO(i,j),obj.jO(i,j)) = v(i,j);
                    end
                end
            end
        end

        function value = blNormGrad(obj,prop)
            q = obj.oGridProp(prop);
            value = (q(:,2)-q(:,1))./(obj.yBL(:,2)-obj.yBL(:,1));
            value = [value (q(:,3:end)-q(:,1:end-2))./(obj.yBL(:,3:end)-obj.yBL(:,1:end-2))];
            value = [value (q(:,end)-q(:,end-1))./(obj.yBL(:,end)-obj.yBL(:,end-1))];

        end



%         function plot_BL_profile(obj,x,prop,ax)
%             if nargin < 4 || isempty(ax)
%                 ax = gca;
%                 disp('Creating axes')
%             end
% 
%             [q, i] = BLprof(obj,x,prop);
%             if string(prop) == "dsdy"
%                 plot(ax, q, obj.yBL(i,2:end))
%             else
%                 plot(ax, q, obj.yBL(i,:))%/obj.yBL(i,end))
%             end
%         end

        function [i, j, blk] = grid_inds_at_BL_max(obj,prop,x)
            io = obj.x2ind(x);
            prop
            prof = obj.BLprof(x, prop);
            [~, jo] = max(prof);
            i = obj.iO(io, jo);
            j = obj.jO(io, jo);
            blk = obj.blkO(io, jo);
        end


        function ind = x2ind(obj,x)
            [~, ind] = min(abs(obj.xSurf-x));
        end

        function value = x2prop(obj, x, prop)
            vals = obj.(prop);
            i = obj.x2ind(x);
            value = vals(i);
        end

        function inlet_profile(obj, prop)
            [prof, y] = putchwise_profile
        end

        function [value, y] = pitchwise_profile(obj, x, prop)

            for ib = 1:obj.NB
                ymin(ib) = min(obj.blk.y{ib}, [], 'all');
                ymax(ib) = max(obj.blk.y{ib}, [], 'all');
            end

            ymin = min(ymin);
            ymax = max(ymax);

            n = 200;
            y = linspace(ymin, ymax, n);
            x(1:n) = x;

            value = obj.unstructured_sample(x, y, prop);

            [xb, yb] = domain_boundary(obj.blk);
            inds = inpolygon(x,y,xb,yb);
            value = value(inds);
            y = y(inds);

        end

        function plot(obj,prop,ax,lims,label,viewarea,rot)
            
            if nargin < 7 || isempty(rot)
                rot = 0;
            end

            R = [cosd(rot) -sind(rot); sind(rot) cosd(rot)];

            if nargin < 3 || isempty(ax)
                ax = gca;
            end
            q = obj.(prop);
            hold on
            for i=1:obj.NB
                x = obj.blk.x{i};
                y = obj.blk.y{i};
                ni = size(x,1);

                coords = [reshape(x, 1, []); reshape(y, 1, [])];
                coords = R' * coords;

                xnow = reshape(coords(1,:), ni, []);
                ynow = reshape(coords(2,:), ni, []);

                pcolor(ax, xnow , ynow, q{i});
            end
            shading('interp')

           if nargin > 6 && ~isempty(viewarea)
                aspect = [(viewarea(2)-viewarea(1)) (viewarea(4)-viewarea(3)) 1];
                pbaspect(aspect)
                axis(viewarea); 
            elseif ~isempty(obj.blk.viewarea)
                pbaspect([(obj.blk.viewarea(2)-obj.blk.viewarea(1)) ...
                    (obj.blk.viewarea(4)-obj.blk.viewarea(3)) 1]);
                axis(obj.blk.viewarea);
            end
            axis equal

            if string(prop) == "schlieren"
                colormap(gray)
                map = colormap;
                map = flip(map,1);
                colormap(map);
                if nargin < 6
                    label = '$|\nabla \rho|/\rho$';
                end
            end

            cb = colorbar;
            if nargin > 3 && ~isempty(lims)
                caxis(lims)
            end
            if nargin > 4 && exist("label","var")
                cb.Label.Interpreter = 'latex';
                cb.Label.String = label;
            end

            set(ax, 'FontSize', 12)
        end

        function blPlot(obj,prop,ax,lims)
            if nargin < 3 || isempty(ax)
                ax = gca;
            end
            q = obj.(prop);
            hold on
            for i=1:obj.NB
                pcolor(ax, obj.blk.x{i}, obj.blk.y{i}, q{i});
            end
            shading('interp')
            axis([-0.1 1.1 0 0.15])
            axis equal
            cb = colorbar;
            if nargin == 5
                caxis(lims)
            end
        end

        function [xsurf, ysurf, n] = getSurfCoordsNorms(obj)
            xo = [];
            yo = [];
            for i=1:length(obj.oblocks)
                nb = obj.oblocks(i);
                xtmp = obj.blk.x{nb}(:,:);
                ytmp = obj.blk.y{nb}(:,:);
                xtmp = flip(xtmp,2);
                ytmp = flip(ytmp,2);
                if obj.oblocks_flip(i) == 1
                    xtmp = flip(xtmp);
                    ytmp = flip(ytmp);
                end
                xo = [xo; xtmp(1:end-1,:)];
                yo = [yo; ytmp(1:end-1,:)];
            end
            xsurf = xo(:,1);
            ysurf = yo(:,1);

            %obj.xSurf = xsurf([obj.iLE:-1:1 end:-1:obj.iTE]);
            %size(obj.xSurf)
            R = [0 -1; 1 0];
            for i = 1:size(xo,1)
                
                if i==1
                    s1 = [(xo(i+1,1)-xo(i,1)); (yo(i+1,1)-yo(i,1))];
                    s2 = [(xo(i,1)-xo(end,1)); (yo(i,1)-yo(end,1))];
                elseif i==size(xo,1)
                    s1 = [(xo(1,1)-xo(i,1)); (yo(1,1)-yo(i,1))];
                    s2 = [(xo(i,1)-xo(i-1,1)); (yo(i,1)-yo(i-1,1))];
                else
                    s1 = [(xo(i+1,1)-xo(i,1)); (yo(i+1,1)-yo(i,1))];
                    s2 = [(xo(i,1)-xo(i-1,1)); (yo(i,1)-yo(i-1,1))];
                end
                n1 = R*s1/norm(s1);
                n2 = R*s2/norm(s2);
                nnow = (0.5*(n1+n2)/norm(0.5*(n1+n2)));
                n(:,i) = nnow;
            end
        end

        function vol = slice2flowBlock(obj, ib)
            vol = volFlowBlock;
            blknow.x = obj.blk.x{ib};
            blknow.y = obj.blk.y{ib};
            vol.blk = blknow;
            vol.ro = obj.ro{ib};
            vol.u = obj.u{ib};
            vol.v = obj.v{ib};
            vol.w = obj.w{ib};
            vol.Et = obj.Et{ib};
            vol.ib = ib;
            vol.blk.nk = 1;
        end

        function newFlow = interpOntoNewGrid(obj, newcase)
            
            newFlow = eval([class(obj) '(newcase.blk, newcase.gas, newcase.bcs)']);
            p = properties(obj);

            xv = [];
            yv = [];

            for ib = 1:obj.NB
                xv = [xv; reshape(obj.blk.x{ib}, [], 1)];
                yv = [yv; reshape(obj.blk.y{ib}, [], 1)];
            end

            for i=1:length(p)
                if iscell(obj.(p{i})) && length(obj.(p{i})) == obj.NB
                    fprintf('Interpolating %s\n', p{i})
                    dv = [];
                    for ib = 1:obj.NB
                        dv = [dv; reshape(obj.(p{i}){ib}, [], 1)];
                    end
                    inds = ~isnan(dv);
                    di = scatteredInterpolant(xv(inds),yv(inds),dv(inds),'linear','boundary');
                    for ib = 1:obj.NB
                        newFlow.(p{i}){ib} = di(newcase.blk.x{ib}, newcase.blk.y{ib});
                    end
                end
            end

        end

        function value = align_tensor_with_surface(obj, prop)

            q = obj.oGridProp(prop);
            value = zeros(size(q));
            
            for i=1:size(q,1)

                n = obj.n(:,i);
                R = [n(2) -n(1)   0; ...
                     n(1)  n(2)   0; ...
                     0     0      1];

                for j=1:size(q,2)
                    value(i,j,:,:) = R*squeeze(q(i,j,:,:))*R';
                end


            end

        end

    end
end