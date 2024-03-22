function slice2kPlot(slice, blk, prop, fpath, lims, label, area, aspect)

    repeats = 1;
    pitch = 0;
    if isfield(blk, "n_pitchwise_repeats")
        repeats = blk.n_pitchwise_repeats;
        pitch = blk.pitch;
    end
    h = figure('Visible','off');
    ax = axes(h);
  
    q = slice.(prop);

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
    pbaspect(aspect)
    axis(area)
    axis off
    if ismember(string(prop),["vortZ","v","w"])
        colormap(redblue)
    end
    cb = colorbar(ax);
    if nargin > 4 && ~isempty(lims)
        caxis(lims);
    end
    if nargin > 5 && ~isempty(label)
        cb.Label.String = label;
    end
    set(ax,'FontSize',12);
    exportgraphics(h, fpath, 'Resolution', 600);
    
    close(h)
    clear
end
