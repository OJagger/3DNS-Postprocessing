function write_case(casename,blk,bcs,gas,solver,nkproc, varargin)

p = inputParser;
addParameter(p, 'casetype', 'both');
parse(p, varargin{:})
nargin
if nargin < 6 || ~exist('nkproc','var')
    nkproc = ceil(solver.nk/solver.npp);
    fprintf('k procs: %d\n', nkproc)
end

if ~isfield(gas,'gamma')
    gas.gamma = gas.gam;
end

NB = length(blk.x);
ncorner = length(blk.corner);   

dir = fullfile(pwd,casename);

if(~exist(dir,'dir'))
mkdir(dir);
end

write_input_files(casename,blk,bcs,gas,solver,'nkproc', nkproc, 'casetype', p.Results.casetype)

% write grid
for ii=1:NB
x = blk.x{ii};
y = blk.y{ii};
[ni,nj]=size(x);
fid = fopen(fullfile(dir,['grid_',num2str(ii),'.txt']),'w');
for j=1:nj
for i=1:ni
fprintf(fid,'%20.16e %20.16e\n',x(i,j),y(i,j));
end
end
fclose(fid);
end


% write blockdims file
fid = fopen(fullfile(dir,'blockdims.txt'),'w');
for ii=1:NB
x = blk.x{ii};
[ni,nj]=size(x);
fprintf(fid,'%d %d %d\n',ni,nj,solver.nk);
end
fclose(fid);
end


