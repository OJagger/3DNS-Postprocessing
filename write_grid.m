function write_grid(path,blk)

% dir = fullfile(pwd,casename);
fprintf('Writing grid files to directory: %s\n',path)

if(~exist(path,'dir'))
mkdir(path);
end

NB = length(blk.x);


% write grid
for ii=1:NB
x = blk.x{ii};
y = blk.y{ii};
[ni,nj]=size(x);
fid = fopen(fullfile(path,['grid_',num2str(ii),'.txt']),'w');
for j=1:nj
for i=1:ni
fprintf(fid,'%20.16e %20.16e\n',x(i,j),y(i,j));
end
end
fclose(fid);
end


% write blockdims file
fid = fopen(fullfile(path,'blockdims.txt'),'w');
for ii=1:NB
x = blk.x{ii};
[ni,nj]=size(x);
fprintf(fid,'%d %d %d\n',ni,nj,blk.nk);
end
fclose(fid);
end