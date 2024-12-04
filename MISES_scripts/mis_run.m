function q = mis_run(directory,run_type)
% MIS_RUN  Run MISES at the design point or a loss loop
%
%   q = MIS_RUN(directory,run_type)
%
%   directory - string designating path to MISES run directory
%   run_type - string to determine run type
%
%   run_type is string
%       'init' - to initialise the grid before running anything with iset
%       'des' - to run a calculation at the design incidence
%       'stall' - to run a loop in the positive incidence direction
%       'loop' - to run a loss loop in both directions

% Default to design point
if exist('run_type','var') == 0
    run_type = 'des';
end

% Record MATLAB directory
mat_dir = pwd;

% Construct Python script directoty
py_dir = [fullfile(mat_dir,'MIS','scripts') '/'];

% Change to MISES directory and set environment variables
cd(directory)
setenv('GFORTRAN_STDIN_UNIT', '5') 
setenv('GFORTRAN_STDOUT_UNIT', '6') 
setenv('GFORTRAN_STDERR_UNIT', '0')
setenv('LD_LIBRARY_PATH','/opt/intel/fce/10.1.015/lib:/home/jvt24/lib:/home/jvt24/lib')

% Include a directory change for windows machines
if isunix == 0
    unix_dir = strrep(strrep(directory,'C:\Users\compressors','~'),'\','/');
    dir_change = ['cd ' unix_dir ' && '];
    py_dir = strrep(strrep(py_dir,'C:\Users\compressors','~'),'\','/');
else
    dir_change = '';
end

% Run MISES functions 
if strcmp(run_type,'init') == 1
    
    % Run the iset command to initialise the grid and solution
    [~,~] = system(['bash --login -c "' dir_change 'python ' py_dir 'mis_run_iset.py"']);
    
    % Return an empty array
    q = [];
    
elseif strcmp(run_type,'des') == 1
    
    % Run the ises command to converge at the design incidence
    [~,~] = system(['bash --login -c "' dir_change 'python ' py_dir 'mis_run_ises.py"']);
    
    % Read the output from the idat file
    q = mis_read_idat('mises',directory);
    
elseif strcmp(run_type,'stall') == 1
    
    % Constants for running loss loops
    loss_lim = 1.25; dAlpha = 1;
    
    % Read the design point loss
    q{1,1} = mis_read_idat('mises',directory);
    loss_des = q{1}.omeg; loss = loss_des;
    alpha_des = atand(q{1}.sinl); alpha = alpha_des;
    
    % Run the positive incidence loop with fixed steps in incidence
    n = 1;
    while loss(end) < loss_des * loss_lim && isempty(q{end}) == 0
    
        % Copy the old idat file to a new name
        copyfile('idat.mises',['idat.mises_' num2str(n)]);
        n = n+1;
    
        % Update the ises input file with new inlet angle
        alpha(n,1) = alpha(end) + dAlpha;
        update_ises(alpha(end));
        
        % Run ises at the current incidence
        [~,~] = system(['bash --login -c "' dir_change 'python ' py_dir filesep 'mis_run_ises.py"']);
        
        % Read the output file and update the loss
        q{end+1,1} = mis_read_idat('mises',directory);
        loss(n,1) = q{end}.omeg;
        
    end

    % Check the loop has crossed the limit
    if loss(end) > loss_des * loss_lim
    
        % Find where the loss crosses the limit
        dl = loss - loss_lim * loss_des;
        j = find(dl(2:end) .* dl(1:end-1) < 0); j = j(end);

        % Interpolate the inlet angle close to failure
        js = j-1:j+1; js(js < 1 | js > length(loss)) = []; 
        alpha_rng = interp1(loss(js),alpha(js),loss_lim * loss_des,'pchip');
        
        % Run MISES again close to the limit
        copyfile('idat.mises',['idat.mises_' num2str(n)]);
        update_ises(alpha_rng);
        [~,~] = system(['bash --login -c "' dir_change 'python ' py_dir filesep 'mis_run_ises.py"']);
        q{end+1,1} = mis_read_idat('mises',directory);
        
    end
    
end

% Change back to the MATLAB home directory
cd(mat_dir)


end


function [] = update_ises(Alpha)
% Update the inlet angle for running loss loops with ises

% Read in the old ises file
A = fileread('ises.mises');

% Update the inlet angle
B = regexprep(A,'\S*(\s\S*\s\S*\s\WMinl)',[sprintf('%6.5f',tand(Alpha)) '$1']);

% Write out the new ises file
fid = fopen('ises.mises','w');
fprintf(fid,B);
fclose(fid);

end

