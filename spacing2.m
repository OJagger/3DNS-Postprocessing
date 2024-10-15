function [r, nj] = spacing2(ywall, y1, dx, target_aspect)
    
    nj = y1/ywall;
    a = dx/ywall;

    while a > target_aspect
        nj = nj-1;
        r = fexpan(y1/ywall, nj);
        a = dx/(ywall*r^(nj-1));
    end
end