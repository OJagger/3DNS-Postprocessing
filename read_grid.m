function [blk, blocks] = read_grid(case_name)
    blk = {};

    tmp = split(case_name, '/');
    if length(tmp) > 1
        case_name = tmp{end};
        base_folder = fullfile('/', tmp{1:end-1});
    else
        base_folder = pwd;
    end

    path = fullfile(base_folder, case_name);
    blocks = readmatrix(fullfile(path, 'blockdims.txt'));
    if size(blocks,1) == 1
        fid = fopen(fullfile(path, 'blockdims.txt'));
        blocks = str2num(char(split(fgetl(fid))))';
        fclose(fid);
    end
    NB = size(blocks,1);
    blk.blockdims = zeros(NB,3);
    for i=1:NB
        [ni,nj,nk] = feval(@(x) x{:}, num2cell(blocks(i,:)));
        blk.blockdims(i,:) = [ni nj nk];
        grid = readmatrix(fullfile(path, ['grid_' num2str(i) '.txt']));
        blk.x{i} = reshape(grid(:,1),[ni,nj]);
        blk.y{i} = reshape(grid(:,2),[ni,nj]);
    end
    blk.nk = nk;
end