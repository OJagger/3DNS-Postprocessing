function [x, y] = domain_boundary(blk)

    ib = 1;

    nb = ib;
    pt = [];

    segments = {};

    for ib=1:length(blk.x)
        pt = get_bnd_patch(blk, ib);

        for p = pt
            switch p
                case "im"
                    segments{end+1} = [reshape(blk.x{ib}(1,:), [], 1) ...
                        reshape(blk.y{ib}(1,:), [], 1)];
                case "ip"
                    segments{end+1} = [reshape(blk.x{ib}(end,:), [], 1) ...
                        reshape(blk.y{ib}(end,:), [], 1)];
                case "jm"
                    segments{end+1} = [reshape(blk.x{ib}(:,1), [], 1) ...
                        reshape(blk.y{ib}(:,1), [], 1)];
                case "jp"
                    segments{end+1} = [reshape(blk.x{ib}(:,end), [], 1) ...
                        reshape(blk.y{ib}(:,end), [], 1)];
            end
        end
    end

    x = segments{1}(:,1);
    y = segments{1}(:,2);
    segments = segments(2:end);

    while x(1) ~= x(end) || y(1) ~= y(end)

        for i=1:length(segments)
            if segments{i}(1,1) == x(end) && segments{i}(1,2) == y(end)
                x = [x; segments{i}(:,1)];
                y = [y; segments{i}(:,2)];
                segments = segments(1:length(segments) ~= i);
                break
            elseif segments{i}(end,1) == x(end) && segments{i}(end,2) == y(end)
                x = [x; segments{i}(end:-1:1,1)];
                y = [y; segments{i}(end:-1:1,2)];
                segments = segments(1:length(segments) ~= i);    
                break
            end
        end

    end

end

function p = get_bnd_patch(blk, ib)

    p = [];
    if blk.next_block{ib}.im == 0
        p = [p, "im"];
    else
        nb = blk.next_block{ib}.im;
        ni = blk.blockdims(nb,1);
        nj = blk.blockdims(nb,2);
        switch blk.next_patch{ib}.im
            case 1
                i=1;
                j=1;
            case 2
                i=ni;
                j=1;
            case 3
                i=1;
                j=1;
            case 4
                i=1;
                j=nj;
        end

        if blk.x{ib}(1,1) ~= blk.x{nb}(i,j) || blk.y{ib}(1,1) ~= blk.y{nb}(i,j)
            p = [p, "im"];
        end
    end

    if blk.next_block{ib}.ip == 0
        p = [p, "ip"];
    else
        nb = blk.next_block{ib}.ip;
        ni = blk.blockdims(nb,1);
        nj = blk.blockdims(nb,2);
        switch blk.next_patch{ib}.ip
            case 1
                i=1;
                j=1;
            case 2
                i=ni;
                j=1;
            case 3
                i=1;
                j=1;
            case 4
                i=1;
                j=nj;
        end

        if blk.x{ib}(end,1) ~= blk.x{nb}(i,j) || blk.y{ib}(end,1) ~= blk.y{nb}(i,j)
            p = [p, "ip"];
        end
    end

    if blk.next_block{ib}.jm == 0
        p = [p, "jm"];
    else
        nb = blk.next_block{ib}.jm;
        ni = blk.blockdims(nb,1);
        nj = blk.blockdims(nb,2);
        switch blk.next_patch{ib}.jm
            case 1
                i=1;
                j=1;
            case 2
                i=ni;
                j=1;
            case 3
                i=1;
                j=1;
            case 4
                i=1;
                j=nj;
        end

        if blk.x{ib}(1,1) ~= blk.x{nb}(i,j) || blk.y{ib}(1,1) ~= blk.y{nb}(i,j)
            p = [p, "jm"];
        end
    end

    if blk.next_block{ib}.jp == 0
        p = [p, "jp"];
    else
        nb = blk.next_block{ib}.jp;
        ni = blk.blockdims(nb,1);
        nj = blk.blockdims(nb,2);
        switch blk.next_patch{ib}.jp
            case 1
                i=1;
                j=1;
            case 2
                i=ni;
                j=1;
            case 3
                i=1;
                j=1;
            case 4
                i=1;
                j=nj;
        end

        if blk.x{ib}(1,end) ~= blk.x{nb}(i,j) || blk.y{ib}(1,end) ~= blk.y{nb}(i,j)
            p = [p, "jp"];
        end
    end


end
