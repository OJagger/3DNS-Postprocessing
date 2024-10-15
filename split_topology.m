function [blk, next_block2, next_patch2, corner2] = split_topology(blk, next_block, next_patch, splits, corner)

    % Function to split 3DNS meshes into smaller blocks. Takes "splits"
    % array as input, with size(NB,2). The ibth row says how many blocks to
    % split block ib into in the i and j directions: eg [2 2] would split
    % the block into 4.

    NB = length(blk.x);
    nbnow = 1;

    % Find new block nos

    for ib=1:NB
        n=prod(splits(ib,:),2)-1;
        tmp = reshape(nbnow:nbnow+n, splits(ib,1), splits(ib,2));
        new_blks{ib}.all = nbnow:nbnow+n;
        new_blks{ib}.im = tmp(1,:);
        new_blks{ib}.ip = tmp(end,:);
        new_blks{ib}.jm = tmp(:,1);
        new_blks{ib}.jp = tmp(:,end);
        nbnow = nbnow+n+1;
    end


    % Check splits are compatible and update next_block, next_patch
    % for ib = 1:NB
    %     for d = dirs
    %         nb = next_block{ib}.(d);
    %         np = next_patch{ib}.(d);
    %         if nb~=0
    %             switch np
    %                 case {1, 2}
    %                     if splits(ib,1) == splits(nb,1)
    %                         for 
    %                 case {3, 4}
    %             end
    %         end
    %     end
    % end


    nd(1:2) = 2;
    nd(3:4) = 1;
    dirs = ["im", "ip", "jm", "jp"];

    % Split blocks and find new patches

    for ib = NB:-1:1

        % if prod(splits(ib,:)) > 1

            xtmp = blk.x{ib};
            ytmp = blk.y{ib};
            [ni, nj] = size(xtmp);
            nbi = splits(ib,1);
            nbj = splits(ib,2);
            nbn = nbi*nbj;
            blk.x = [blk.x(1:ib) cell(1, nbn-1) blk.x(ib+1:end)];
            blk.y = [blk.y(1:ib) cell(1, nbn-1) blk.y(ib+1:end)];
            % next_patch2 = [next_patch2(1:ib) cell(1, nbn-1) next_patch2(ib+1:end)];
            % next_block2 = [next_block2(1:ib) cell(1, nbn-1) next_block2(ib+1:end)];
            % blk.x{ib+nbn:end+nbn-1} = blk.x{ib+1:end};
            % blk.y{ib+nbn:end+nbn-1}  = blk.y{ib+1:end};
    
            ii = round(linspace(1, ni, nbi+1));
            jj = round(linspace(1, nj, nbj+1));

            ibnow = ib;
            ibnew = new_blks{ib}.all(1);

            for j=1:nbj
                for i=1:nbi
                    ibnow;
                    blk.x{ibnow} = xtmp(ii(i):ii(i+1),jj(j):jj(j+1));
                    blk.y{ibnow} = ytmp(ii(i):ii(i+1),jj(j):jj(j+1));
                    ibnow = ibnow+1;

                    if i == 1
                        nb = next_block{ib}.im;
                        np = next_patch{ib}.im;

                        if nb == 0
                            next_block2{ibnew}.im = 0;
                        else
                            if splits(ib,2) ~= splits(nb,nd(np))
                                error(sprintf('Block %d and %d have incompatible splits in i direction', [ib, nb]));
                            end
                            next_block2{ibnew}.im = new_blks{nb}.(dirs(np))(j);
                        end
                        next_patch2{ibnew}.im = np;

                    else
                        next_block2{ibnew}.im = ibnew-1;
                        next_patch2{ibnew}.im = 2;
                    end
                    
                    if i == nbi
                        nb = next_block{ib}.ip;
                        np = next_patch{ib}.ip;
                        


                        if nb == 0
                            next_block2{ibnew}.ip = 0;
                        else
                            if splits(ib,2) ~= splits(nb,nd(np))
                                error(sprintf('Block %d and %d have incompatible splits in i direction', [ib, nb]));
                            end
                            next_block2{ibnew}.ip = new_blks{nb}.(dirs(np))(j);
                        end
                        next_patch2{ibnew}.ip = np;

                    else
                        next_block2{ibnew}.ip = ibnew+1;
                        next_patch2{ibnew}.ip = 1;
                    end

                    if j == 1
                        nb = next_block{ib}.jm;
                        np = next_patch{ib}.jm;

                        
                        if nb == 0
                            next_block2{ibnew}.jm = 0;
                        else
                            if splits(ib,1) ~= splits(nb,nd(np))
                                error(sprintf('Block %d and %d have incompatible splits in j direction', [ib, nb]));
                            end
                            next_block2{ibnew}.jm = new_blks{nb}.(dirs(np))(i);
                        end
                        next_patch2{ibnew}.jm = np;

                    else
                        next_block2{ibnew}.jm = ibnew-nbi;
                        next_patch2{ibnew}.jm = 4;
                    end

                    if j == nbj
                        nb = next_block{ib}.jp;
                        np = next_patch{ib}.jp;
                        
                        
                        if nb == 0
                            next_block2{ibnew}.jp = 0;
                        else
                            if splits(ib,1) ~= splits(nb,nd(np))
                                error(sprintf('Block %d and %d have incompatible splits in j direction', [ib, nb]));
                            end
                            next_block2{ibnew}.jp = new_blks{nb}.(dirs(np))(i);
                        end
                        next_patch2{ibnew}.jp = np;

                    else
                        next_block2{ibnew}.jp = ibnew+nbi;
                        next_patch2{ibnew}.jp = 3;
                    end
                    
                    ibnew = ibnew+1;
                end
            end

        % end


    end

    % Now update corner struct with new block nos

    for ic = 1:length(corner)
        corner2{ic}.Nb = corner{ic}.Nb;
        for k=1:corner{ic}.Nb
            ib = corner{ic}.block{k};
            ni = blk.blockdims(ib,1);
            nj = blk.blockdims(ib,2);
            if corner{ic}.i{k} == 1 && corner{ic}.j{k} == 1
                corner2{ic}.block{k} = new_blks{ib}.im(1);
                corner2{ic}.i{k} = 1;
                corner2{ic}.j{k} = 1;
            elseif corner{ic}.i{k} > 1 && corner{ic}.j{k} == 1
                corner2{ic}.block{k} = new_blks{ib}.ip(1);
                corner2{ic}.i{k} = size(blk.x{corner2{ic}.block{k}},1);
                corner2{ic}.j{k} = 1;
            elseif corner{ic}.i{k} == 1 && corner{ic}.j{k} > 1
                corner2{ic}.block{k} = new_blks{ib}.im(end);
                corner2{ic}.i{k} = 1;
                corner2{ic}.j{k} = size(blk.x{corner2{ic}.block{k}},2);
            elseif corner{ic}.i{k} > 1 && corner{ic}.j{k} >1
                corner2{ic}.block{k} = new_blks{ib}.ip(end);
                corner2{ic}.i{k} = size(blk.x{corner2{ic}.block{k}},1);
                corner2{ic}.j{k} = size(blk.x{corner2{ic}.block{k}},2);
            end
        end
    end

    for ib = 1:length(blk.x)
        blk.blockdims(ib,1:2) = size(blk.x{ib});
        blk.procdims(ib,:) = round(blk.blockdims(ib,:)/blk.npp);
        blk.blockdims(ib,3) = blk.nk;
    end

    oblocks = blk.oblocks;
    oblocks_flip = blk.oblocks_flip;
    blk.oblocks = [];
    blk.oblocks_flip = [];
    for i = 1:length(oblocks)
        ib = oblocks(i);
        new_ibs = new_blks{ib}.all;
        blk.oblocks = [blk.oblocks new_ibs];
        blk.oblocks_flip = [blk.oblocks_flip oblocks_flip(i)*ones(size(new_ibs))];
    end
end