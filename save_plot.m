
function save_plot(dir, fname, save, fmt, draft)

    if nargin < 3 || isempty(save)
        save = true;
    end
    if nargin < 4
        fmt = 'png';
    end
    if nargin < 5
        draft = false;
    end
    fmts = string(fmt);
    if save
        f = gcf;
        for fmt = fmts
            switch fmt
                case "png"

                    % exportgraphics(f, fullfile(dir, [fname '.png']), 'Resolution',128)
                    exportgraphics(f, fullfile(dir, [fname '.png']), 'Resolution',1000)
                case "eps"
                    exportgraphics(f, fullfile(dir, [fname '.eps']))
                case "svg"
                    set(f, 'Color', 'None');
                    figure(f)
                    set(gca, 'Color', 'None');
                    plot2svg(fullfile(dir, [fname '.svg']), f, 'png')
                case "pdf"
                    exportgraphics(f, fullfile(dir, [fname '.pdf']))
                case "fig"
            end
        end
%        savefig(f, fullfile(dir, [fname '.fig']))
    end

end
