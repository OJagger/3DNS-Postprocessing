function blkNodes = writeFluentCas(path, blk, boundaries, iWrite)

    if nargin < 4
        iWrite = true;
    end

    meshPath = strrep(path,'.cas','.msh');
    writeFluentMesh(meshPath, blk, boundaries, iWrite);
    
    jouPath = fullfile(fileparts(path),'mesh2cas.jou');

    f = fopen(jouPath,'w');
    fprintf(f, ...
        ['file/read-case %s\n' ...
        'file/cff-files no\n' ...
        'file/write-case %s\n' ...
        '/exit'], meshPath, path);
    fclose(f);

    system(['fluent -r23.1.0 2ddp -g t1 -i' jouPath]);

end
