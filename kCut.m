classdef kCut < flowSlice
    %KCUT Superclass for instantaneous and mean flow on a k boundary

    properties
        xSurf;
        yBL;
        xO;
        yO;
        iO;
        jO;
        blkO;
        oblocks;
        oblocks_flip;
        iLE;
        iTE;
        n;
        ssurf;          % Surface distance fron LE
        vortZ;          % Z vorticity
    end

    methods
        function obj = kCut(blk, gas)
            obj@flowSlice(blk,gas)
            disp('Constucting kCut')
            if nargin > 0
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
                        ni = size(blk.x{obj.oblocks(i)},1);
                        nj = size(blk.x{obj.oblocks(i)},2);
                        xtmp = blk.x{obj.oblocks(i)}(:,:);
                        ytmp = blk.y{obj.oblocks(i)}(:,:);
                        itmp = 1:ni;
                        itmp = repmat(itmp',[1 nj]);
                        jtmp = 1:nj;
                        jtmp = repmat(flip(jtmp), [ni 1]);
                        blktmp = obj.oblocks(i)*ones(ni,nj);
                        xtmp = flip(xtmp,2);
                        ytmp = flip(ytmp,2);
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
                        end

                    end
                    xsurf = xo(:,1);
                    [~, obj.iLE] = min(xsurf);
                    [~, obj.iTE] = max(xsurf);
                    obj.xSurf = xsurf(obj.iLE:obj.iTE);
                    obj.xO = xo(obj.iLE:obj.iTE,:);
                    obj.yO = yo(obj.iLE:obj.iTE,:);
                    obj.iO = io(obj.iLE:obj.iTE,:);
                    obj.jO = jo(obj.iLE:obj.iTE,:);
                    obj.blkO = blko(obj.iLE:obj.iTE,:);
                    %obj.xSurf = xsurf([obj.iLE:-1:1 end:-1:obj.iTE]);
                    %size(obj.xSurf)
                    R = [0 -1; 1 0];

                    obj.yBL = zeros(obj.iTE-obj.iLE+1,size(xo,2));
                    obj.ssurf = zeros(1,size(obj.yBL,1));
                    for i = obj.iLE:obj.iTE
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
                        obj.n(:,i+1-obj.iLE) = nnow;
                        for j=1:size(xo,2)
                            dx = xo(i,j) - xo(i,1);
                            dy = yo(i,j) - yo(i,1);
                            obj.yBL(i+1-obj.iLE,j) = dot(nnow, [dx;dy]);
                        end
                        if i>obj.iLE
                            dx = xo(i,1) - xo(i-1,1);
                            dy = yo(i,1) - yo(i-1,1);
                            ds = sqrt(dx^2 + dy^2);
                            obj.ssurf(i+1-obj.iLE) = obj.ssurf(i-obj.iLE) + ds;
                        end
                    end
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


        function [profile, i] = BLprof(obj, x, prop)
            %BLPROF get a profile of prop across the BL at specified x
            [~, i] = min(abs(obj.xSurf-x));
            if any(strcmp(["dsdy" "U"],prop))
                profile = obj.(prop)(i,:);
                size(profile)
            else
                propfield = oGridProp(obj,prop);
                profile = propfield(i,:);
            end
        end

        function blfield = oGridProp(obj, prop)
            %OGRIDPROP Construct array of property in o grid for one surface
            %of blade
            
            if any(strcmp(["U","dsdy"], prop))
                blfield = obj.(prop);
            else
                blfield = [];
                propnow = obj.(prop);
                for i=1:length(obj.oblocks)
                    clear temp
                    temp = propnow{obj.oblocks(i)};
                    temp = flip(temp,2);
                    %size(temp)
                    if obj.oblocks_flip(i) == 1
                        temp = flip(temp);
                    end
                    blfield = [blfield; temp(1:end-1,:)];
                end
                blfield = blfield(obj.iLE:obj.iTE,:);
                %size(blfield,1)
                %blfield = blfield([obj.iLE:-1:1 end:-1:obj.iTE],:);
                %size(blfield,1)
            end
        end

        function value = blNormGrad(obj,prop)
            q = obj.oGridProp(prop);
            value = (q(:,2)-q(:,1))./(obj.yBL(:,2)-obj.yBL(:,1));
            value = [value (q(:,3:end)-q(:,1:end-2))./(obj.yBL(:,3:end)-obj.yBL(:,1:end-2))];
            value = [value (q(:,end)-q(:,end-1))./(obj.yBL(:,end)-obj.yBL(:,end-1))];

        end

        function plot_BL_profile(obj,x,prop,ax)
            if nargin < 4 || isempty(ax)
                ax = gca;
                disp('Creating axes')
            end

            [q, i] = BLprof(obj,x,prop);
            if string(prop) == "dsdy"
                plot(ax, q, obj.yBL(i,2:end))
            else
                plot(ax, q, obj.yBL(i,:)/obj.yBL(i,end))
            end
        end

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


        

        function plot(obj,prop,ax,lims,label)
            
            if nargin < 3 || isempty(ax)
                ax = gca;
            end
            q = obj.(prop);
            hold on
            for i=1:obj.NB
                pcolor(ax, obj.blk.x{i}, obj.blk.y{i}, q{i});
            end
            shading('interp')
            
            if ~isempty(obj.blk.viewarea)
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
                xtmp = obj.blk.x{obj.oblocks(i)}(:,:);
                ytmp = obj.blk.y{obj.oblocks(i)}(:,:);
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
    end
end