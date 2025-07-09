function [i, r] = in_domain(blk, x, y, pitch)

    if nargin < 4
        pitch = blk.pitch;
    end

    [xb, yb] = domain_boundary(blk);

    if nargout == 1
        i = inpolygon(x, y, xb, yb);

    else
        for r=-4:4
            i = inpolygon(x, y-r*pitch, xb, yb);
            if i
                break
            end
        end
        if ~i
            r = [];
        end
    end

end