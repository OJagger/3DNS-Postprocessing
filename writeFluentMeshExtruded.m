function blkNodes = writeFluentMeshExtruded(path, blk, boundaries, span, nk, iWrite)

    % Write 3DNS style prismatic mesh as Fluent .msh file in ASCII format
    % path: path to write file (ending in .msh)
    % blk: struct containing mesh. Must contain fields:
    %
    %   x: Cell array of x coord arrays (size (ni, nj)), one per block
    %   y: Cell array of y coord arrays (size (ni, nj)), one per block
    %
    %   next_block: Cell array of structs containing fields with keys im, 
    %       ip, jm, jp. The field values define the block connectivity. 
    %       Eg. next_block{3}.ip = 5 means the i=ni face of block 3 
    %       boarders block 5. Set to zero on boundaries.
    %
    %   next_patch: Cell array of structs containing fields with keys im, 
    %       ip, jm, jp. The values define which face of the neighboring
    %       block is connected:
    %       1: im (i=1) face
    %       2: ip (i=ni) face
    %       3: jm (j=1) face
    %       4: jp (j=nj) face
    %       Eg. next_patch{3}.ip = 3 means that the i=ni face of block 3
    %       boarders the j=1 face of block 5 (following previous example)
    %       On boundaries, the values define the boundary types in 3DNS:
    %       1 for inlets
    %       2 for outlets
    %       3 for walls
    %       however these boundary values are not needed here
    %
    %   boundries: Cell array of structs defining the boundaries. Each has
    %       fields:
    %       label: string label for boundary - written into file (for
    %       debugging only)
    %       blocks: list of blocks contained in boundary
    %       patches: list defining which face of each block makes up
    %       boundary. 1 for im, 2 for ip etc.
    %       type. Boundary type (Fluent spec):
    %       3: Wall
    %       4: Pressure inlet
    %       5: Pressure outlet
    %       8: Periodic shadow
    %       10: Velocity inlet
    %       12: Periodic
    %       see romeo.univ-reims.fr/documents/fluent/tgrid/ug/appb.pdf
    %       for complete list
    %
    % span: spanwise extent (extrusion distance)
    % nk: numper of grid points in span (No desired cells + 1)
    % iWrite: set to false to disable file writeout for debugging
    
    


    next_block = blk.next_block;
    next_patch = blk.next_patch;

    if nargin < 7
        iWrite = true;
    end
    
    NB = length(blk.x);
    
    im_assigned(1:NB) = false;
    ip_assigned(1:NB) = false;
    jm_assigned(1:NB) = false;
    jp_assigned(1:NB) = false;
    
    nodes = [];
    faces = [];
    cells = [];
    nnodes = 1;
    ncells = 1;
    nfaces = 1;

    zv = linspace(0,span,nk);

    for ib = 1:NB
        [ni, nj] = size(blk.x{ib});
        blk.nodes{ib} = zeros(ni, nj, nk);
        blk.ifaces{ib} = zeros(ni, nj-1, nk-1);
        blk.jfaces{ib} = zeros(ni-1, nj, nk-1);
        blk.kfaces{ib} = zeros(ni-1, nj-1, nk);
    end
    
    for ib = 1:NB
    
        is = 1:size(blk.x{ib},1);
        js = 1:size(blk.x{ib},2);
    
        if im_assigned(ib)
            blk.nodes{ib}(1,:,:) = getBoarderNodes(ib, 'im');
            is = is(2:end);
        end
        if ip_assigned(ib)
            blk.nodes{ib}(end,:,:) = getBoarderNodes(ib, 'ip');
            is = is(1:end-1);
        end
        if jm_assigned(ib)
            blk.nodes{ib}(:,1,:) = getBoarderNodes(ib, 'jm');
            js = js(2:end);
        end
        if jp_assigned(ib)
            blk.nodes{ib}(:,end,:) = getBoarderNodes(ib, 'jp');
            js = js(1:end-1);
        end
    
        for k=1:nk
            for i=is
                for j=js
                    blk.nodes{ib}(i,j,k) = nnodes;
                    nodes(nnodes).x = blk.x{ib}(i,j);
                    nodes(nnodes).y = blk.y{ib}(i,j);
                    nodes(nnodes).z = zv(k);
                    nnodes = nnodes+1;
                end
            end
        end

        setAssigned(ib, 'im');
        setAssigned(ib, 'ip');
        setAssigned(ib, 'jm');
        setAssigned(ib, 'jp');
    
    end

    % Set cell Nos

    for ib = 1:NB
        for k=1:nk-1
            for i = 1:size(blk.x{ib},1)-1
                for j = 1:size(blk.x{ib},2)-1
                    blk.cells{ib}(i,j,k) = ncells;
                    ncells = ncells+1;
                end
            end
        end
    end

    % Set face Nos
    
    im_assigned(1:NB) = false;
    ip_assigned(1:NB) = false;
    jm_assigned(1:NB) = false;
    jp_assigned(1:NB) = false;
    boundary_zone(1:NB,1:6) = 0;

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
                for k=1:nk-1
                    for j = 1:size(blk.x{ib},2)-1
                        blk.ifaces{ib}(1,j,k) = nfaces;
                        faces(nfaces).bc = boundaries{b}.type;
                        nfaces = nfaces+1;
                    end
                end
            elseif p == 2
                for k=1:nk-1
                    for j = 1:size(blk.x{ib},2)-1
                        blk.ifaces{ib}(end,j,k) = nfaces;
                        faces(nfaces).bc = boundaries{b}.type;
                        nfaces = nfaces+1;
                    end
                end
            elseif p == 3
                for k=1:nk-1
                    for i = 1:size(blk.x{ib},1)-1
                        blk.jfaces{ib}(i,1,k) = nfaces;
                        faces(nfaces).bc = boundaries{b}.type;
                        nfaces = nfaces+1;
                    end
                end
            else
                for k=1:nk-1
                    for i = 1:size(blk.x{ib},1)-1
                        blk.jfaces{ib}(i,end,k) = nfaces;
                        faces(nfaces).bc = boundaries{b}.type;
                        nfaces = nfaces+1;
                    end
                end
            end
        end
        boundaries{b}.nEnd = nfaces-1;
    end

    % zm peroiodic
    boundaries{end+1}.nStart = nfaces;
    for ib = 1:NB
        boundary_zone(ib,5) = length(boundaries);
        next_block{ib}.km = ib;
        next_block{ib}.kp = ib;
        next_patch{ib}.km = 6;
        next_patch{ib}.kp = 5;
        for i = 1:size(blk.x{ib},1)-1
            for j = 1:size(blk.x{ib},2)-1
    %             faces(nfaces).n1 = blk.nodes{ib}(end,j);
    %             faces(nfaces).n2 = blk.nodes{ib}(end,j+1);
                blk.kfaces{ib}(i,j,1) = nfaces;
                faces(nfaces).bc = 12;
                nfaces = nfaces+1;
            end
        end
    end
    boundaries{end}.label = "km Periodic";
    boundaries{end}.type = 12;
    boundaries{end}.blocks = 1:NB;
    boundaries{end}.patches(1:NB) = 5;
    boundaries{end}.nEnd = nfaces-1;
    

    % zp peroiodic
    boundaries{end+1}.nStart = nfaces;
    for ib = 1:NB
        boundary_zone(ib,6) = length(boundaries);
        for i = 1:size(blk.x{ib},1)-1
            for j = 1:size(blk.x{ib},2)-1
    %             faces(nfaces).n1 = blk.nodes{ib}(end,j);
    %             faces(nfaces).n2 = blk.nodes{ib}(end,j+1);
                blk.kfaces{ib}(i,j,end) = nfaces;
                faces(nfaces).bc = 8;
                nfaces = nfaces+1;
            end
        end
    end
    boundaries{end}.label = "kp Periodic";
    boundaries{end}.type = 8;
    boundaries{end}.blocks = 1:NB;
    boundaries{end}.patches(1:NB) = 6;
    boundaries{end}.nEnd = nfaces-1;


    for ib = 1:NB
        fprintf('Calculating block %d/%d\n',ib,NB)
        [ni, nj] = size(blk.x{ib});
        if ib == 1
            ni;
        end

        % Set remaining boarder faces
        if boundary_zone(ib,1) == 0
            if ~im_assigned(ib)
                for k=1:nk-1
                    for j = 1:nj-1
                        blk.ifaces{ib}(1,j,k) = nfaces;
                        nfaces = nfaces+1;
                    end
                end
                setAssigned(ib, 'im')
            else
                blk.ifaces{ib}(1,:,:) = getBoarderFaces(ib,'im');
            end
        end

        if boundary_zone(ib,2) == 0
            if ~ip_assigned(ib)
                for k=1:nk-1
                    for j = 1:nj-1
                        blk.ifaces{ib}(end,j,k) = nfaces;
                        nfaces = nfaces+1;
                    end
                end
                setAssigned(ib, 'ip');
            else
                blk.ifaces{ib}(end,:,:) = getBoarderFaces(ib,'ip');
            end
        end

        if boundary_zone(ib,3) == 0
            if ~jm_assigned(ib)
                for k=1:nk-1
                    for i = 1:ni-1
                        blk.jfaces{ib}(i,1,k) = nfaces;
                        nfaces = nfaces+1;
                    end
                end
                setAssigned(ib, 'jm');
            else
                blk.jfaces{ib}(:,1,:) = getBoarderFaces(ib,'jm');
            end
        end

        if boundary_zone(ib,4) == 0
            if ~jp_assigned(ib)
                for k=1:nk-1
                    for i = 1:ni-1
                        blk.jfaces{ib}(i,end,k) = nfaces;
                        nfaces = nfaces+1;
                    end
                end
                setAssigned(ib, 'jp');
            else
                blk.jfaces{ib}(:,end,:) = getBoarderFaces(ib,'jp');
            end
        end

        % Set interior face Nos

        % i faces
        for k=1:nk-1
            for i=2:ni-1
                for j=1:nj-1
                    blk.ifaces{ib}(i,j,k) = nfaces;
                    nfaces = nfaces+1;
                end
            end
        end
        % j faces
        for k=1:nk-1
            for j=2:nj-1
                for i=1:ni-1
                    blk.jfaces{ib}(i,j,k) = nfaces;
                    nfaces = nfaces+1;
                end
            end
        end
        % k faces
        for k=2:nk-1
            for j=1:nj-1
                for i=1:ni-1
                    blk.kfaces{ib}(i,j,k) = nfaces;
                    nfaces = nfaces+1;
                end
            end
        end

        % Now get face nodes and cells
        
        % i faces
        r1 = [blk.x{ib}(2,1) - blk.x{ib}(1,1); ...
            blk.y{ib}(2,1) - blk.y{ib}(1,1); 0];
        r2 = [blk.x{ib}(1,2) - blk.x{ib}(1,1); ...
            blk.y{ib}(1,2) - blk.y{ib}(1,1); 0];
        r3 = [0; 0; zv(2)-zv(1)];
        swap = false;
        if(dot(cross(r1,r2),r3)) < 0
            swap = true;
        end
        for i=1:ni
            mcells = zeros(1,nj-1,nk-1);
            pcells = zeros(1,nj-1,nk-1);
            if i==1
                mcells = getBoarderCells(ib, 'im');
                pcells = blk.cells{ib}(i,:,:);
            elseif i==ni
                mcells = blk.cells{ib}(i-1,:,:);
                pcells = getBoarderCells(ib, 'ip');
            else
                mcells = blk.cells{ib}(i-1,:,:);
                pcells = blk.cells{ib}(i,:,:);
            end

            if swap
                tmp = mcells;
                mcells = pcells;
                pcells = tmp;
            end

            for k=1:nk-1
                for j=1:nj-1
                    n = blk.ifaces{ib}(i,j,k);
                    if (i==1) && (j==1) && (k==2)
                        disp('')
                    end
                    faces(n).n1 = blk.nodes{ib}(i,j,k);
                    faces(n).n2 = blk.nodes{ib}(i,j+1,k);
                    faces(n).n3 = blk.nodes{ib}(i,j+1,k+1);
                    faces(n).n4 = blk.nodes{ib}(i,j,k+1);
                    faces(n).cr = pcells(1,j,k);
                    faces(n).cl = mcells(1,j,k);
                    if faces(n).bc == 0
                        faces(n).bc = 2;
                    end
                end
            end
        end

        % j faces
%         r1 = [blk.x{ib}(2,1) - blk.x{ib}(1,1); ...
%             blk.y{ib}(2,1) - blk.y{ib}(1,1); 0];
%         r2 = [blk.x{ib}(1,2) - blk.x{ib}(1,1); ...
%             blk.y{ib}(1,2) - blk.y{ib}(1,1); 0];
%         r3 = [0; 0; zv(2)-zv(1)];
%         swap = false;
        if(dot(cross(r1,r2),r3)) < 0
            swap = true;
        end
        for j=1:nj
            mcells = zeros(ni-1,1,nk-1);
            pcells = zeros(ni-1,1,nk-1);
            if j==1
                mcells = blk.cells{ib}(:,j,:);
                pcells = getBoarderCells(ib, 'jm');
            elseif j==nj
                mcells = getBoarderCells(ib, 'jp');
                pcells = blk.cells{ib}(:,j-1,:);
            else
                mcells = blk.cells{ib}(:,j,:);
                pcells = blk.cells{ib}(:,j-1,:);
            end

            if swap
                tmp = mcells;
                mcells = pcells;
                pcells = tmp;
            end

            for k=1:nk-1
                for i=1:ni-1
                    n = blk.jfaces{ib}(i,j,k);
                    if n == 0
                        n
                    end
                    faces(n).n1 = blk.nodes{ib}(i,j,k);
                    faces(n).n2 = blk.nodes{ib}(i+1,j,k);
                    faces(n).n3 = blk.nodes{ib}(i+1,j,k+1);
                    faces(n).n4 = blk.nodes{ib}(i,j,k+1);
                    faces(n).cr = pcells(i,1,k);
                    faces(n).cl = mcells(i,1,k);
                    if faces(n).bc == 0
                        faces(n).bc = 2;
                    end
                end
            end
        end

        % k faces
        for k=1:nk
            mcells = zeros(ni-1,nj-1,1);
            pcells = zeros(ni-1,nj-1,1);
            if k==1
                pcells = blk.cells{ib}(:,:,k);
                %pcells = getBoarderCells(ib, 'jm');
            elseif k==nk
                %mcells = getBoarderCells(ib, 'jp');
                mcells = blk.cells{ib}(:,:,k-1);
            else
                pcells = blk.cells{ib}(:,:,k);
                mcells = blk.cells{ib}(:,:,k-1);
            end

            if swap
                tmp = mcells;
                mcells = pcells;
                pcells = tmp;
            end

            for i=1:ni-1
                for j=1:nj-1
                    n = blk.kfaces{ib}(i,j,k);
                    if n == 0
                        n
                    end
                    faces(n).n1 = blk.nodes{ib}(i,j,k);
                    faces(n).n2 = blk.nodes{ib}(i+1,j,k);
                    faces(n).n3 = blk.nodes{ib}(i+1,j+1,k);
                    faces(n).n4 = blk.nodes{ib}(i,j+1,k);
                    faces(n).cr = pcells(i,j,1);
                    faces(n).cl = mcells(i,j,1);
                    if faces(n).bc == 0
                        faces(n).bc = 2;
                    end
                    if k==nk
                        faces(n).sf = blk.kfaces{ib}(i,j,1);
                    end
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
                case 5
                    dr = 'km';
                case 6
                    dr = 'kp';
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
                        per_faces = blk.ifaces{ib}(1,:,:);
                    case 2
                        dr = 'ip';
                        per_faces = blk.ifaces{ib}(end,:,:);
                    case 3
                        dr = 'jm';
                        per_faces = blk.jfaces{ib}(:,1,:);
                    case 4
                        dr = 'jp';
                        per_faces = blk.jfaces{ib}(:,end,:);
                    case 5
                        dr = 'km';
                        per_faces = blk.kfaces{ib}(:,:,1);
                    case 6
                        dr = 'kp';
                        per_faces = blk.kfaces{ib}(:,:,end);

                end

                if p == 5
                    shadow_faces = blk.kfaces{ib}(:,:,end);
                elseif p == 6
                    shadow_faces = blk.kfaces{ib}(:,:,1);
                else
                    shadow_faces = getBoarderFaces(ib, dr);
                end

                per_faces = squeeze(per_faces);
                shadow_faces = squeeze(shadow_faces);
                
                for i = 1:size(per_faces,1)
                    for j = 1:size(per_faces,2)
                        pairs(npairs).face = per_faces(i,j);
                        pairs(npairs).shadow_face = shadow_faces(i,j);
                        npairs = npairs+1;
                    end
                end
            end
            boundaries{b}.nPairsEnd = npairs-1;

        end
    end

%     fprintf('nnodes: %d, should be %d\n', length(nodes), ni*nj*nk)
%     fprintf('ncells: %d, should be %d\n', ncells-1, ((ni-1)*(nj-1)*(nk-1)))
%     fprintf('nfaces: %d, should be %d\n', length(faces), ((ni-1)*nj*(nk-1) + ni*(nj-1)*(nk-1) + (ni-1)*(nj-1)*nk))

    facesStr = cell(length(faces),1);
    nodesStr = cell(length(nodes),1);
    shadowStr = cell(1+boundaries{end}.nEnd-boundaries{end}.nStart,1);

    fprintf('Creating ASCII face information to write\n')
    parfor n=1:length(faces)
        facesStr{n} = sprintf('%X %X %X %X %X %X',[faces(n).n1 faces(n).n2, faces(n).n3 faces(n).n4, ...
                faces(n).cr faces(n).cl]);
    end

    fprintf('Creating ASCII node information to write\n')
    parfor n=1:length(nodes)
        nodesStr{n} = sprintf('%10.8f %10.8f %10.8f', [nodes(n).x nodes(n).y nodes(n).z]);
    end

%     fprintf('Creating ASCII shadow zone information\n')
%     for i=1:length(shadowStr)
%         n = i-1+boundaries{end}.nStart;
%         shadowStr{i} = sprintf('%X %X',[faces(n).sf n]);
%     end

    blkNodes = blk.nodes;

    if iWrite
        
        path;
        fid = fopen(path,'w');
    
        fprintf(fid,'(0 "Grid:")\n');
        fprintf(fid,'\n');
    
        fprintf(fid,'(0 "Dimensions:")\n');
        fprintf(fid,'(2 3)\n');
        fprintf(fid,'\n');
        
        fprintf(fid,'(0 "Cells:")\n');
        fprintf(fid,'(12 (0 %X %X 0))\n',[1 (ncells-1)]);
        fprintf(fid, '\n');
    
        fprintf(fid,'(0 "Faces:")\n');
        fprintf(fid,'(13 (0 %X %X 0))\n',[1 length(faces)]);
        fprintf(fid, '\n');
    
        fprintf(fid,'(0 "Nodes:")\n');
        fprintf(fid,'(10 (0 %X %X 0 3))\n',[1 length(nodes)]);
        fprintf(fid, '\n');
        
        fprintf(fid,'(12 (b 1 %X 1 4))\n',ncells-1);
        fprintf(fid, '\n');


        nZones = 2;

        if boundaries{end}.nEnd < length(faces)
            fprintf('Writing interior faces\n')
            fprintf(fid,'(0 "Interior:")\n');
            fprintf(fid,'(13 (2 %X %X %X 4)(\n', [boundaries{end}.nEnd+1 length(faces) 2]);
            for n=(boundaries{end}.nEnd+1):length(faces)
                fprintf(fid,'%X %X %X %X %X %X\n',[faces(n).n1 faces(n).n2, faces(n).n3 faces(n).n4, ...
                    faces(n).cr faces(n).cl]);
            end
    %         fprintf(fid, '%s\n', facesStr{boundaries{end}.nEnd+1:end});
            fprintf(fid, '))\n');
            fprintf(fid, '\n');
        else
            nZones = 1;
        end
    
        for ib = 1:length(boundaries)
            b = boundaries{ib};
            nZones = nZones+1;
            boundaries{ib}.nZone = nZones;
            fprintf(fid,'(0 "%s:")\n', b.label);
            fprintf(fid,'(13 (%d %X %X %X 4)(\n', ...
                [nZones b.nStart b.nEnd b.type]);
            for n=b.nStart:b.nEnd
                if faces(n).cr == 0
                    fprintf(fid,'%X %X %X %X %X %X\n',[faces(n).n1 faces(n).n4 faces(n).n3 faces(n).n2 faces(n).cl faces(n).cr]);
                else
                    fprintf(fid,'%X %X %X %X %X %X\n',[faces(n).n1 faces(n).n2 faces(n).n3 faces(n).n4 faces(n).cr faces(n).cl]);
                end
            end
            fprintf(fid, '))\n');
            fprintf(fid, '\n');
        end

        fprintf('Writing shadow zone\n')
        fprintf(fid,'(0 "Shadow zone pairs:")\n');

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
    
        fprintf('Writing nodes\n')
        fprintf(fid,'(10 (1 %X %X 1 3)\n(\n',[1 length(nodes)]);
%         for i=1:length(nodes)
%             fprintf(fid,'%10.8f %10.8f %10.8f\n', [nodes(i).x nodes(i).y nodes(i).z]);
%         end
        fprintf(fid, '%s\n', nodesStr{:});
        fprintf(fid, '))\n');
        fprintf(fid, '\n');
        
        fclose(fid);
    end

    function nodesnow = getBoarderNodes(ib, dr)
        nblk = next_block{ib}.(dr);
        nptch = next_patch{ib}.(dr);
        switch nptch
            case 1
            nodesnow = blk.nodes{nblk}(1,:,:);
            if ismember(dr,{'jm','jp'})
                nodesnow = permute(nodesnow, [2 1 3]);
            end
            case 2
                nodesnow = blk.nodes{nblk}(end,:,:);
            if ismember(dr,{'jm','jp'})
                nodesnow = permute(nodesnow, [2 1 3]);
            end
            case 3
                nodesnow = blk.nodes{nblk}(:,1,:);
            if ismember(dr,{'im','ip'})
                nodesnow = permute(nodesnow, [2 1 3]);
            end
            case 4
                nodesnow = blk.nodes{nblk}(:,end,:);
            if ismember(dr,{'im','ip'})
                nodesnow = permute(nodesnow, [2 1 3]);
            end
        end
    end

    function facesnow = getBoarderFaces(ib, dr)
        nblk = next_block{ib}.(dr);
        nptch = next_patch{ib}.(dr);
            switch nptch
                case 1
                    facesnow = blk.ifaces{nblk}(1,:,:);
                    if ismember(dr,{'jm','jp'})
                        facesnow = permute(facesnow, [2 1 3]);
                    end
                case 2
                    facesnow = blk.ifaces{nblk}(end,:,:);
                    if ismember(dr,{'jm','jp'})
                        facesnow = permute(facesnow, [2 1 3]);
                    end
                case 3
                    facesnow = blk.jfaces{nblk}(:,1,:);
                    if ismember(dr,{'im','ip'})
                        facesnow = permute(facesnow, [2 1 3]);
                    end
                case 4
                    facesnow = blk.jfaces{nblk}(:,end,:);
                    if ismember(dr,{'im','ip'})
                        facesnow = permute(facesnow, [2 1 3]);
                    end
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
                    cellsnow(1,1:njb-1,1:nk-1) = 0;
                case {'jm','jp'}
                    cellsnow(1:nib-1,1,1:nk-1) = 0;
            end
        else
            nptch = next_patch{ib}.(dr);
            switch nptch
                case 1
                    cellsnow = blk.cells{nblk}(1,:,:);
                    if ismember(dr,{'jm','jp'})
                        cellsnow = permute(cellsnow, [2 1 3]);
                    end
                case 2
                    cellsnow = blk.cells{nblk}(end,:,:);
                    if ismember(dr,{'jm','jp'})
                        cellsnow = permute(cellsnow, [2 1 3]);
                    end
                case 3
                    cellsnow = blk.cells{nblk}(:,1,:);
                    if ismember(dr,{'im','ip'})
                        cellsnow = permute(cellsnow, [2 1 3]);
                    end
                case 4
                    cellsnow = blk.cells{nblk}(:,end,:);
                    if ismember(dr,{'im','ip'})
                        cellsnow = permute(cellsnow, [2 1 3]);
                    end
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
