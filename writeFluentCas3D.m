function writeFluentCas3D(path, blk, boundaries)

    % Write 3DNS style prismatic mesh, extruded 1 cell in the span, and
    % convert to .cas
    % path: path to write cas file (ending with .cas)
    % blk: struct containing mesh information
    % boundaries: struct defining boundary patches

    nk = 2;
    span = 1.0; 

    meshPath = strrep(path,'.cas','.msh');
    writeFluentMeshExtruded(meshPath, blk, boundaries, span, nk);
    
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
