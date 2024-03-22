function blkNodes = writeFluentMesh(path, blk, boundaries, iWrite)

    if nargin < 5
        iWrite = true;
    end

    next_block = blk.next_block;
    next_patch = blk.next_patch;
    
    NB = length(blk.x);
    
    im_assigned(1:NB) = false;
    ip_assigned(1:NB) = false;
    jm_assigned(1:NB) = false;
    jp_assigned(1:NB) = false;

    % Set boundary zones to exclude from boarder node and face sharing
%     for b = 1:length(boundaries)
%         for n = 1:length(boundaries{b}.blocks)
%             ib = boundaries{b}.blocks(n);
%             p = boundaries{b}.patches(n);
%             boundary_zone(ib,p) = true;
%         end
%     end
    
    nodes = [];
    faces = [];
    cells = [];
    nnodes = 1;
    ncells = 1;
    nfaces = 1;

    for ib = 1:NB
        [ni, nj] = size(blk.x{ib});
        blk.nodes{ib} = zeros(ni, nj);
        blk.ifaces{ib} = zeros(ni, nj-1);
        blk.jfaces{ib} = zeros(ni-1, nj);
    end
    
    for ib = 1:NB
    
        is = 1:size(blk.x{ib},1);
        js = 1:size(blk.x{ib},2);
    
        if im_assigned(ib)
            blk.nodes{ib}(1,:) = getBoarderNodes(next_block{ib}.im, next_patch{ib}.im);
            is = is(2:end);
        end
        if ip_assigned(ib)
            blk.nodes{ib}(end,:) = getBoarderNodes(next_block{ib}.ip, next_patch{ib}.ip);
            is = is(1:end-1);
        end
        if jm_assigned(ib)
            blk.nodes{ib}(:,1) = getBoarderNodes(next_block{ib}.jm, next_patch{ib}.jm);
            js = js(2:end);
        end
        if jp_assigned(ib)
            blk.nodes{ib}(:,end) = getBoarderNodes(next_block{ib}.jp, next_patch{ib}.jp);
            js = js(1:end-1);
        end
    
        for i=is
            for j=js
                blk.nodes{ib}(i,j) = nnodes;
                nodes(nnodes).x = blk.x{ib}(i,j);
                nodes(nnodes).y = blk.y{ib}(i,j);
                nnodes = nnodes+1;
            end
        end

        setAssigned(ib, 'im');
        setAssigned(ib, 'ip');
        setAssigned(ib, 'jm');
        setAssigned(ib, 'jp');
    
    end

    % Set cell nos

    for ib = 1:NB
        for i = 1:size(blk.x{ib},1)-1
            for j = 1:size(blk.x{ib},2)-1
                blk.cells{ib}(i,j) = ncells;
                ncells = ncells+1;
            end
        end
    end

    % Set face nos
    
    im_assigned(1:NB) = false;
    ip_assigned(1:NB) = false;
    jm_assigned(1:NB) = false;
    jp_assigned(1:NB) = false;
    boundary_zone(1:NB,1:4) = 0;

    for b = 1:length(boundaries)
        boundaries{b}.nStart = nfaces;
        for n = 1:length(boundaries{b}.blocks)
            ib = boundaries{b}.blocks(n);
            p = boundaries{b}.patches(n);
            boundary_zone(ib,p) = b;

            % BCs:
            % 3: Wall
            % 4: Pressure inlet
            % 5: Pressure outlet
            % 12: Periodic
            % 8: Periodic shadow

            if p == 1
                for j = 1:size(blk.x{ib},2)-1
                    blk.ifaces{ib}(1,j) = nfaces;
                    faces(nfaces).bc = boundaries{b}.type;
                    nfaces = nfaces+1;
                end
            elseif p == 2
                for j = 1:size(blk.x{ib},2)-1
                    blk.ifaces{ib}(end,j) = nfaces;
                    faces(nfaces).bc = boundaries{b}.type;
                    nfaces = nfaces+1;
                end
            elseif p == 3
                for i = 1:size(blk.x{ib},1)-1
                    blk.jfaces{ib}(i,1) = nfaces;
                    faces(nfaces).bc = boundaries{b}.type;
                    nfaces = nfaces+1;
                end
            else
                for i = 1:size(blk.x{ib},1)-1
                    blk.jfaces{ib}(i,end) = nfaces;
                    faces(nfaces).bc = boundaries{b}.type;
                    nfaces = nfaces+1;
                end
            end
        end
        boundaries{b}.nEnd = nfaces-1;
    end

    for ib = 1:NB
        fprintf('Calculating block %d/%d\n',ib,NB)
        [ni, nj] = size(blk.x{ib});
        if ib == 1
            ni;
        end

        % Set remaining boarder faces
        if boundary_zone(ib,1) == 0
            if ~im_assigned(ib)
                for j = 1:nj-1
                    blk.ifaces{ib}(1,j) = nfaces;
                    nfaces = nfaces+1;
                end
                setAssigned(ib, 'im')
            else
                blk.ifaces{ib}(1,:) = getBoarderFaces(ib,'im');
            end
        end

        if boundary_zone(ib,2) == 0
            if ~ip_assigned(ib)
                for j = 1:nj-1
                    blk.ifaces{ib}(end,j) = nfaces;
                    nfaces = nfaces+1;
                end
                setAssigned(ib, 'ip');
            else
                blk.ifaces{ib}(end,:) = getBoarderFaces(ib,'ip');
            end
        end

        if boundary_zone(ib,3) == 0
            if ~jm_assigned(ib)
                for i = 1:ni-1
                    blk.jfaces{ib}(i,1) = nfaces;
                    nfaces = nfaces+1;
                end
                setAssigned(ib, 'jm');
            else
                blk.jfaces{ib}(:,1) = getBoarderFaces(ib,'jm');
            end
        end

        if boundary_zone(ib,4) == 0
            if ~jp_assigned(ib)
                for i = 1:ni-1
                    blk.jfaces{ib}(i,end) = nfaces;
                    nfaces = nfaces+1;
                end
                setAssigned(ib, 'jp');
            else
                blk.jfaces{ib}(:,end) = getBoarderFaces(ib,'jp');
            end
        end

        % Set interior face nos

        % i faces
        for i=2:ni-1
            for j=1:nj-1
                blk.ifaces{ib}(i,j) = nfaces;
                nfaces = nfaces+1;
            end
        end
        % j faces
        for j=2:nj-1
            for i=1:ni-1
                blk.jfaces{ib}(i,j) = nfaces;
                nfaces = nfaces+1;
            end
        end

        % Now get face nodes and cells

        r1 = [blk.x{ib}(2,1) - blk.x{ib}(1,1); ...
            blk.y{ib}(2,1) - blk.y{ib}(1,1); 0];
        r2 = [blk.x{ib}(1,2) - blk.x{ib}(1,1); ...
            blk.y{ib}(1,2) - blk.y{ib}(1,1); 0];
        k = [0; 0; 1];
        swap = false;
        if dot(k,cross(r1,r2)) < 0
            swap = true;
        end
        
        % i faces

        for i=1:ni
            rcells = zeros(1,nj-1);
            lcells = zeros(1,nj-1);
            if i==1
                rcells = getBoarderCells(ib, 'im');
                lcells = blk.cells{ib}(i,:);
            elseif i==ni
                rcells = blk.cells{ib}(i-1,:);
                lcells = getBoarderCells(ib, 'ip');
            else
                rcells = blk.cells{ib}(i-1,:);
                lcells = blk.cells{ib}(i,:);
            end

            if swap
                tmp = rcells;
                rcells = lcells;
                lcells = tmp;
            end

            for j=1:nj-1
                n = blk.ifaces{ib}(i,j);
                faces(n).n1 = blk.nodes{ib}(i,j);
                faces(n).n2 = blk.nodes{ib}(i,j+1);
                faces(n).cr = rcells(j);
                faces(n).cl = lcells(j);
                if faces(n).bc == 0
                    faces(n).bc = 2;
                end
            end
        end

        % j faces

        for j=1:nj
            rcells = zeros(1,nj-1);
            lcells = zeros(1,nj-1);
            if j==1
                rcells = blk.cells{ib}(:,j);
                lcells = getBoarderCells(ib, 'jm');
            elseif j==nj
                rcells = getBoarderCells(ib, 'jp');
                lcells = blk.cells{ib}(:,j-1);
            else
                rcells = blk.cells{ib}(:,j);
                lcells = blk.cells{ib}(:,j-1);
            end

            if swap
                tmp = rcells;
                rcells = lcells;
                lcells = tmp;
            end

            for i=1:ni-1
                n = blk.jfaces{ib}(i,j);
                faces(n).n1 = blk.nodes{ib}(i,j);
                faces(n).n2 = blk.nodes{ib}(i+1,j);
                faces(n).cr = rcells(i);
                faces(n).cl = lcells(i);
                if faces(n).bc == 0
                    faces(n).bc = 2;
                end
            end
        end
    end

    % Match periodic faces with shadow zone faces

    npairs = 1;
    for b = 1:length(boundaries)
        if boundaries{b}.type == 12
            ib = boundaries{b}.blocks(1);
            p = boundaries{b}.patches(1);
            switch p
                case 1
                    dr = 'im';
                case 2
                    dr = 'ip';
                case 3
                    dr = 'jm';
                case 4
                    dr = 'jp';
            end
            
            shadow_boundary = boundary_zone(next_block{ib}.(dr),...
                next_patch{ib}.(dr));
            boundaries{b}.shadow_boundary = shadow_boundary;
            boundaries{shadow_boundary}.type = 8;
            boundaries{b}.nPairsStart = npairs;

            for n = 1:length(boundaries{b}.blocks)
                ib = boundaries{b}.blocks(n);
                p = boundaries{b}.patches(n);
                switch p
                    case 1
                        dr = 'im';
                        per_faces = blk.ifaces{ib}(1,:);
                    case 2
                        dr = 'ip';
                        per_faces = blk.ifaces{ib}(end,:);
                    case 3
                        dr = 'jm';
                        per_faces = blk.jfaces{ib}(:,1);
                    case 4
                        dr = 'jp';
                        per_faces = blk.jfaces{ib}(:,end);
                end
                shadow_faces = getBoarderFaces(ib, dr);
                for i = 1:length(per_faces)
                    pairs(npairs).face = per_faces(i);
                    pairs(npairs).shadow_face = shadow_faces(i);
                    npairs = npairs+1;
                end
            end
            boundaries{b}.nPairsEnd = npairs-1;

        end
    end

    fprintf('nnodes: %d, should be %d\n', length(nodes), ni*nj)
    fprintf('ncells: %d, should be %d\n', ncells-1, ((ni-1)*(nj-1)))
    fprintf('nfaces: %d, should be %d\n', length(faces), ((ni-1)*nj + ni*(nj-1)))

    blkNodes = blk.nodes;

    if iWrite
        
        path
        fid = fopen(path,'w');
    
        fprintf(fid,'(0 "Grid:")\n');
        fprintf(fid,'\n');
    
        fprintf(fid,'(0 "Dimensions:")\n');
        fprintf(fid,'(2 2)\n');
        fprintf(fid,'\n');
        
        fprintf(fid,'(0 "Cells:")\n');
        fprintf(fid,'(12 (0 %X %X 0))\n',[1 (ncells-1)]);
        fprintf(fid, '\n');
    
        fprintf(fid,'(0 "Faces:")\n');
        fprintf(fid,'(13 (0 %X %X 0))\n',[1 length(faces)]);
        fprintf(fid, '\n');
    
        fprintf(fid,'(0 "Nodes:")\n');
        fprintf(fid,'(10 (0 %X %X 0 2))\n',[1 length(nodes)]);
        fprintf(fid, '\n');
        
        fprintf(fid,'(12 (9 1 %X 1 3))\n',ncells-1);
        fprintf(fid, '\n');
    
        fprintf(fid,'(0 "Interior:")\n');
        fprintf(fid,'(13 (2 %X %X %X 2)(\n', [boundaries{end}.nEnd+1 length(faces) 2]);
        for n=(boundaries{end}.nEnd+1):length(faces)
            fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2 faces(n).cr faces(n).cl]);
        end
        fprintf(fid, '))\n');
        fprintf(fid, '\n');

        nZones = 2;
        for ib = 1:length(boundaries)
            b = boundaries{ib};
            nZones = nZones+1;
            boundaries{ib}.nZone = nZones;
            fprintf(fid,'(0 "%s:")\n', b.label);
            fprintf(fid,'(13 (%d %X %X %X 2)(\n', ...
                [nZones b.nStart b.nEnd b.type]);
            for n=b.nStart:b.nEnd
                if faces(n).cr == 0
                    fprintf(fid,'%X %X %X %X\n',[faces(n).n2 faces(n).n1 faces(n).cl faces(n).cr]);
                else
                    fprintf(fid,'%X %X %X %X\n',[faces(n).n1 faces(n).n2 faces(n).cr faces(n).cl]);
                end
            end
            fprintf(fid, '))\n');
            fprintf(fid, '\n');
        end

        for ib = 1:length(boundaries)
            b = boundaries{ib};
            if b.type == 12  % Periodic zone
                periodicZone = b.nZone;
                shadowZone = boundaries{b.shadow_boundary}.nZone;
                fprintf(fid,'(0 "%s shadow zone:")\n', b.label);
                fprintf(fid,'(18 (%X %X %X %X)(\n', ...
                [b.nPairsStart b.nPairsEnd ...
                periodicZone shadowZone]);
                for n=b.nPairsStart:b.nPairsEnd
                    fprintf(fid,'%X %X\n',[pairs(n).face pairs(n).shadow_face]);
                end
                fprintf(fid, '))\n');
                fprintf(fid, '\n');
            end
        end
        
    
        fprintf(fid,'(10 (1 %X %X 1 2)\n(\n',[1 length(nodes)]);
        for i=1:length(nodes)
            fprintf(fid,'%10.8f %10.8f\n', [nodes(i).x nodes(i).y]);
        end
        fprintf(fid, '))\n');
        fprintf(fid, '\n');
        
        fclose(fid);
    end

    function nodesnow = getBoarderNodes(nblk, nptch)
        switch nptch
            case 1
                nodesnow = blk.nodes{nblk}(1,:);
            case 2
                nodesnow = blk.nodes{nblk}(end,:);
            case 3
                nodesnow = blk.nodes{nblk}(:,1);
            case 4
                nodesnow = blk.nodes{nblk}(:,end);
        end
    end

    function facesnow = getBoarderFaces(ib, dr)
        nblk = next_block{ib}.(dr);
        nptch = next_patch{ib}.(dr);
            switch nptch
                case 1
                    facesnow = blk.ifaces{nblk}(1,:);
                case 2
                    facesnow = blk.ifaces{nblk}(end,:);
                case 3
                    facesnow = blk.jfaces{nblk}(:,1);
                case 4
                    facesnow = blk.jfaces{nblk}(:,end);
            end
        end

    function cellsnow = getBoarderCells(ib, dr)
        nblk = next_block{ib}.(dr);
        if nblk > NB
            nblk = 0;
        end
        if (nblk == 0) || is_periodic(ib, dr)
            [nib, njb] = size(blk.x{ib});
            switch dr
                case {'im','ip'}
                    cellsnow(1:njb-1) = 0;
                case {'jm','jp'}
                    cellsnow(1:nib-1) = 0;
            end
        else
            nptch = next_patch{ib}.(dr);
            switch nptch
                case 1
                    cellsnow = blk.cells{nblk}(1,:);
                case 2
                    cellsnow = blk.cells{nblk}(end,:);
                case 3
                    cellsnow = blk.cells{nblk}(:,1);
                case 4
                    cellsnow = blk.cells{nblk}(:,end);
            end
        end
    end

    function setAssigned(ib, dr)
        nblk = next_block{ib}.(dr);
        if nblk > NB
            nblk = 0;
        end
        if nblk ~= 0 && ~is_periodic(ib, dr)
            nptch = next_patch{ib}.(dr);
            switch nptch
                case 1
                    im_assigned(nblk) = true;
                case 2
                    ip_assigned(nblk) = true;
                case 3
                    jm_assigned(nblk) = true;
                case 4
                    jp_assigned(nblk) = true;
            end
        end
    end

    function periodic = is_periodic(ib, dr)
        nblk = next_block{ib}.(dr);
        nptch = next_patch{ib}.(dr);
        if nblk ~= 0
            switch dr
                case 'im'
                    x1 = blk.x{ib}(1,1);
                    y1 = blk.y{ib}(1,1);
                case 'ip'
                    x1 = blk.x{ib}(end,1);
                    y1 = blk.y{ib}(end,1);
                case 'jm'
                    x1 = blk.x{ib}(1,1);
                    y1 = blk.y{ib}(1,1);
                case 'jp'
                    x1 = blk.x{ib}(1,end);
                    y1 = blk.y{ib}(1,end);
            end
            switch nptch
                case 1
                    x2 = blk.x{nblk}(1,1);
                    y2 = blk.y{nblk}(1,1);
                case 2
                    x2 = blk.x{nblk}(end,1);
                    y2 = blk.y{nblk}(end,1);
                case 3
                    x2 = blk.x{nblk}(1,1);
                    y2 = blk.y{nblk}(1,1);
                case 4
                    x2 = blk.x{nblk}(1,end);
                    y2 = blk.y{nblk}(1,end);
            end
        end

        if (nblk ~= 0) && ((x1 ~= x2) || (y1 ~= y2))
            periodic = true;
        else
            periodic = false;
        end
    end


end
