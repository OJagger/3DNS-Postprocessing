function h = slice2schlierenVortOverlay(slice, blk, fpath, lims, area, aspect)
    if length(lims) == 1
        lims = abs(lims)*[-1 1];
    end
    h = figure('Visible','off');
    ax1 = axes(h);
    q1 = slice.schlieren;
    pcolor(ax1, blk.x{1}, blk.y{1}, q1{1});
    caxis([0 100]);
    shading interp
    axis equal
    pbaspect(ax1,aspect);
    map = colormap(gray);
    map = flip(map,1);
    axis(ax1,area);
    axis off

    ax2 = axes(h);
    om = abs(slice.vortZ{1});
    a =2000;
    b = 2000;
    mask = 0.5*(1+tanh((om-a)/b));

    q2 = slice.vortZ;
    pcolor(ax2, blk.x{1}, blk.y{1}, q2{1});
    caxis(lims);
    shading interp
    axis equal
    pbaspect(ax2,aspect);
    axis(ax2,area);
    axis off
    alpha(ax2,mask)
    colormap(ax1,map);
    colormap(ax2, redblue);

    set(ax2, 'Position', get(ax1, 'Position'));

    if ~strcmp(fpath,'')
        exportgraphics(h, fpath, 'Resolution', 600);
    end
end
