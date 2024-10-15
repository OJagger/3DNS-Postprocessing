function slice2kPlot(slice, blk, prop, fpath, lims, label, area, aspect)

    repeats = 1;
    pitch = 0;
    if isfield(blk, "n_pitchwise_repeats")
        repeats = blk.n_pitchwise_repeats;
        pitch = blk.pitch;
    end
    h = figure('Visible','off');
    ax = axes(h);
  
    if strcmp(string(prop), "overlay")
        q = slice.schlieren;
        overlay = true;
    else
        q = slice.(prop);
        overlay = false;
    end

    hold on
    offset = 0;
    if repeats > 2
        offset = -blk.pitch;
    end
    for ir = 1:repeats
    for i=1:slice.NB
        pcolor(ax, blk.x{i}, blk.y{i}+offset+(ir-1)*pitch, q{i});
    end
    end
    shading('interp')
    axis equal
    pbaspect(ax, aspect)
    axis(ax, area)
    axis off
    if overlay
        clim([0 100]);
        map = colormap(gray);
        map = flip(map,1);
        colormap(map);
    else
        if ismember(string(prop),["vortZ","v","w"])
            colormap(redblue)
        end
        cb = colorbar(ax);
        if nargin > 4 && ~isempty(lims)
            clim(lims);
        end
        if nargin > 5 && ~isempty(label)
            cb.Label.String = label;
        end
    end
    set(ax,'FontSize',12);


    if overlay
        ax2 = axes(h);
        q2 = slice.vortZ;

        a =2000;
        b = 2000;

        hold on
        for i=1:slice.NB
            om = abs(q2{i});
            mask = 0.5*(1+tanh((om-a)/b));
            for ir = 1:repeats
                s = pcolor(ax2, blk.x{i}, blk.y{i}+offset+(ir-1)*pitch, q2{i});
                alpha(s, mask);
            end
        end

        shading('interp')
        axis equal
        area = axis(ax);
        aspect = [(area(2)-area(1)) (area(4)-area(3)) 1];
        pbaspect(aspect);
        axis off
        colormap(ax2, redblue);
        clim(lims);
        set(ax2, 'Position', get(ax, 'Position'));

    end

    exportgraphics(h, fpath, 'Resolution', 600);
    
    close(h)
    clear
end
