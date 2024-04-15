function blkNodes = writeFluentCas3D(path, blk, boundaries, iWrite)

    if nargin < 4
        iWrite = true;
    end

    nk = 2;
    span = 1.0; % (blk.x{1}(end,1)-blk.x{1}(1,1))*(nk-1)/(size(blk.x{1},1)-1);

    meshPath = strrep(path,'.cas','.msh');
    writeFluentMeshExtruded(meshPath, blk, boundaries, span, nk, iWrite);
    
    jouPath = fullfile(fileparts(path),'mesh2cas.jou');

    f = fopen(jouPath,'w');
    fprintf(f, ...
        ['file/read-case %s\n' ...
        'file/cff-files no\n' ...
        'file/write-case %s\n' ...
        '/exit'], meshPath, path);
    fclose(f);

    system(['fluent -r23.1.0 3ddp -g t1 -i' jouPath]);

end
