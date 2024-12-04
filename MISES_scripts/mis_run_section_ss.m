function [p,h,job] = mis_run_section_ss(directory,h,col,plot_stuff,is_turb,run_single,job)
% MIS_RUN_SECTION  Run MISES on a given section geometry
%
%   [p,h = MIS_RUN_SECTION(directory,c,f,h,col,plot_stuff,is_turb)
%
%   directory - string of directory to run in
%   c - input section data structure or coordinates
%   f - input boundary conditions
%   h - optional figure handle
%   col - optional RGB colour vector
%   plot_stuff - 0 or 1 for showing working
%   is_turb - optional transition inputs
%   p - output data structure
%
%   c is either data structure or cell array:
%       if data structure then it must be a full definition to run in
%       bl_construct_section
%       if cell array then cell 1 is xrrt coords, cell 2 is blade count

%% Prepare defaults
if exist('col','var') == 0 || isempty(col) == 1
    col = [0 0 0];
end
if exist('plot_stuff','var') == 0
    plot_stuff = 1;
end
if exist('is_turb','var') == 0
    is_turb = 1;
end

if exist('run_single','var') == 0
    run_single = 0;
end

%% determine target inlet and exit Mach number triangles
if exist('job','var') == 1
    if isfield(job,'f') == 0
        [f]=mises_output_finput(job);
    else
        f=job.f;
        f.Re = 1*10^6;
        f.M = job.Mrel_inlet;
        f.AVDR = job.f.Ax;
    end
end


if f.Alpha < 0; f.Alpha = -f.Alpha; end;
if f.Alpha2 < 0; f.Alpha2 = -f.Alpha2; end;

%% create directory if it does not exist
job.rjm_directory_temp = job.rjm_directory;
directory = strrep(directory,'TURBOSTREAM','MISES');
% reset job.rjm_directory
job.rjm_directory = directory;

% Make directory if required
if exist(directory,'dir') == 0
    mkdir(directory);
end

% Delete old output files from directory
if isempty(dir([directory 'polar*'])) == 0; delete([directory 'polar*']); end;
if isempty(dir([directory 'output*'])) == 0; delete([directory 'output*']); end;
if isempty(dir([directory 'idat*'])) == 0; delete([directory 'idat*']); end;
if isempty(dir([directory 'spec*'])) == 0; delete([directory 'spec*']); end;

%% create section and parameterise section for running on MISES and write it out
write_file = 1;
[mt,c] = mises_create_mt(job,write_file);

%% Write section file
save([directory 'section.mat'],'c','f');

%% Write stream file
fid = fopen([directory 'stream.mises'],'w');

if isfield(job,'Axt')==0
    
        
    % apply radial contraction from 25% chord upstream to 25% chord downstream
    fprintf(fid,'%i\n',0);
    fprintf(fid,'%6.5f %6.5f %6.5f\n',[min(mt(:,1))-0.5 c.r_le/c.tchord 1.00000]);
    fprintf(fid,'%6.5f %6.5f %6.5f\n',[0.0-0.25*c.m_chord c.r_le/c.tchord 1.00000]);
    fprintf(fid,'%6.5f %6.5f %6.5f\n',[0.0-0.25*c.m_chord c.r_le/c.tchord 1.00000]);
    fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1))+0.25*c.m_chord c.r_te/c.tchord f.AVDR]);
    fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1))+0.25*c.m_chord c.r_te/c.tchord f.AVDR]);
    fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1))+0.5 c.r_te/c.tchord f.AVDR]);
    fclose(fid);
    
else

    
    % changed radial contraction up to throat and back again
    fprintf(fid,'%i\n',0);
    fprintf(fid,'%6.5f %6.5f %6.5f\n',[min(mt(:,1))-0.5 c.r_le/c.tchord 1.00000]);
    fprintf(fid,'%6.5f %6.5f %6.5f\n',[0.0-0.25*c.m_chord c.r_le/c.tchord 1.00000]);
    fprintf(fid,'%6.5f %6.5f %6.5f\n',[0.0-0.25*c.m_chord c.r_le/c.tchord 1.00000]);
    fprintf(fid,'%6.5f %6.5f %6.5f\n',[c.xs_sb*max(mt(:,1)) c.r_le/c.tchord job.Axt]);
    fprintf(fid,'%6.5f %6.5f %6.5f\n',[c.xs_sb*max(mt(:,1)) c.r_le/c.tchord job.Axt]);
    fprintf(fid,'%6.5f %6.5f %6.5f\n',[c.xs_sb*max(mt(:,1)) c.r_le/c.tchord job.Axt]);
    % fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1))+0.40*max(mt(:,1)) c.r_te/c.tchord 1]);
    fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1))+0.25*c.m_chord c.r_te/c.tchord f.AVDR]);
    fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1))+0.25*c.m_chord c.r_te/c.tchord f.AVDR]);
    fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1))+0.5 c.r_te/c.tchord f.AVDR]);
    fclose(fid);
    
end

%% Write ises file
fid = fopen([directory 'ises.mises'],'w');

fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
% supersonic
fprintf(fid,'%s\n',' 15 4 3 17 6 !Global constraints');

fprintf(fid,'%6.5f %6.5f %6.5f %6.5f %6.5f ',[f.M f.Pin/f.Poin tand(f.Alpha) -0.66*c.m_chord f.Vtao1]);
fprintf(fid,'%s\n','|Minl p1/p01 Sinl Xinl v1/ao1');
fprintf(fid,'%6.5f %6.5f %6.5f %6.5f  ', [f.M2 f.Pout/f.Poin tand(f.Alpha2) 1.33*c.m_chord]);
fprintf(fid,'%s\n','|Mout p2/p01 Sout Xout');
fprintf(fid,'%s\n','0.00000     0.00000     1.39433     0.00000  | MFR  HWRATin GAMin DIFACTOR');
fprintf(fid,'%i ',round(f.Re));
fprintf(fid,'%s\n', '-4.00000 0.00000  0.26133 |Re -Turb');
if is_turb == 1
    fprintf(fid,'%s\n','0.001 0.001  |Xtr1 Xtr2');
else
    fprintf(fid,'%s\n','1.10 1.10  |Xtr1 Xtr2');
end
fprintf(fid,'%s\n','4     0.65000  -4.00000  0.00000   | ISMOM  MCRIT  MUCON  PLOSSIN');
fprintf(fid,'%s\n','0.00000  0.00000   | BVR1in  BVR2in');
fclose(fid);


%% create mesh
n_type_max = 7;
n_type = 1; % n=100
mis_gridpar(directory,n_type);

% Change to MISES directory and set environment variables
cd(directory)
setenv('GFORTRAN_STDIN_UNIT', '5')
setenv('GFORTRAN_STDOUT_UNIT', '6')
setenv('GFORTRAN_STDERR_UNIT', '0')
setenv('LD_LIBRARY_PATH','/opt/intel/fce/10.1.015/lib:/home/dl467/lib:/home/dl467/lib')

%% Run MISES and change directory back

% supersonic
[~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_ss.py"','-echo');

p.manual = 0;
cd(directory)

% Read in flow file and re-run if not converged
if exist([directory 'polarx.mises'],'file') ~= 0
    [Polarx, Ises] = mis_read_polarx('mises',directory);
    if isstruct(Polarx) ~= 1
        disp('*******Run again*********');
        while isstruct(Polarx) ~= 1
            n_type = n_type + 1;
            
            if n_type<=n_type_max
                disp(['****Create new gridpar ' num2str(n_type) '****']);
                mis_gridpar(directory,n_type);
                
                % supersonic
                [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_ss.py"','-echo');
                
                p.manual = 0;
                [Polarx, Ises] = mis_read_polarx('mises',directory);
            else
                disp('*******Have to run manually*********');
                while isstruct(Polarx) ~= 1
                    if f.M<1
                        [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_manual.py"','-echo');
                    else
                        p.manual = 1;
                        Polarx = [];
                        Polarx.p = {};
                        break
                    end
                    [Polarx, Ises] = mis_read_polarx('mises',directory);
                end
            end
        end
        n_type = 1;
    end
else
    disp('File Not Found')
    disp('*******Run again*********');
    Polarx = [];
    while isstruct(Polarx) ~= 1
        n_type = n_type + 1;
        if n_type<=n_type_max
            disp(['****Create new gridpar ' num2str(n_type) '****']);
            mis_gridpar(directory,n_type);
            
            % supersonic
            [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_ss.py"','-echo');
            p.manual = 0;
            
            [Polarx, Ises] = mis_read_polarx('mises',directory);
        else
            disp('*******Have to run manually*********');
            while isstruct(Polarx) ~= 1
                if f.M<1
                    [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_manual.py"','-echo');
                else
                    p.manual = 1;
                    Polarx = [];
                    Polarx.p = {};
                    
                    break
                    
                end
                [Polarx, Ises] = mis_read_polarx('mises',directory);
            end
        end
    end
    n_type = 1;
    
end

%% determine inlet blade metal angle that will give inlet relative flow angle and exit Mach number constraint

[Polarx, Ises] = mis_read_polarx('mises',directory);
binl_target = Ises.binl;
if p.manual == 1
    disp('*****Cannot run initial iteration*****')
    binl = binl_target;
else
    binl = Polarx.binl;
end

% iteration number
ni = 1;
n_psi_max = 10;
n_type_max = 7;

while abs(binl_target-binl)>0.15 && ni < n_psi_max && isstruct(Polarx) == 1 && p.manual == 0
    
    [~,~]=system(['rm ' directory 'polarx.mises'],'-echo');
    
    dinlet_cam = binl_target-binl;
    
    if abs(dinlet_cam) > 1
        dinlet_cam = 1*sign(dinlet_cam);
    end
    if isfield(job,'chi_le')==1
        chi_le = job.chi_le + dinlet_cam;
        job.chi_le_0 = job.chi_le;
        binl_0 = binl;
        job.chi_le=chi_le;
    else % legay
        inlet_cam = job.inlet_cam + dinlet_cam;
        job.inlet_cam_0 = job.inlet_cam;
        binl_0 = binl;
        job.inlet_cam=inlet_cam;
    end
    
    [mt,c] = mises_create_mt(job,1);
    % Write section file
    save([directory 'section.mat'],'c','f');
    % supersonic
    [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_ss.py"','-echo');
    
    
    ni = ni+1;
    
    % Read in flow file and re-run if not converged
    if exist([directory 'polarx.mises'],'file') ~= 0
        [Polarx, Ises] = mis_read_polarx('mises',directory);
        if isstruct(Polarx) ~= 1
            disp('*******Run again*********');
            while isstruct(Polarx) ~= 1 && n_type < n_type_max + 1
                n_type = n_type + 1;
                disp(['****Create new gridpar ' num2str(n_type) '****']);
                mis_gridpar(directory,n_type);
                
                % supersonic
                [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_ss.py"','-echo');
                [Polarx, Ises] = mis_read_polarx('mises',directory);
                
                p.manual = 0;
            end
            if n_type > n_type_max
                disp('*******Have to run manually*********');
                while isstruct(Polarx) ~= 1
                    p.manual = 1;
                    Polarx = []; Polarx.p = {};
                    break
                    
                end
            end
            n_type = 1; mis_gridpar(directory,n_type);
        end
    else
        disp('File Not Found')
        disp('*******Run again*********');
        Polarx = [];
        while isstruct(Polarx) ~= 1 && n_type < n_type_max + 1
            
            n_type = n_type + 1;
            disp(['****Create new gridpar ' num2str(n_type) '****']);
            mis_gridpar(directory,n_type);
            
            % supersonic
            [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_ss.py"','-echo');
            [Polarx, Ises] = mis_read_polarx('mises',directory);
            
        end
        p.manual = 0;
        if n_type> n_type_max
            disp('*******Have to run manually*********');
            while isstruct(Polarx) ~= 1
                p.manual = 1;
                Polarx = []; Polarx.p = {};
                break
            end
        end
        n_type = 1; mis_gridpar(directory,n_type);
        
    end
    
    if exist([directory 'polarx.mises'],'file') ~= 0 && p.manual == 0
        [Polarx, Ises] = mis_read_polarx('mises',directory);
        binl_target = Ises.binl;
        binl = Polarx.binl;
        display(['Difference in inlet flow angle to target: ' num2str(binl-binl_target)])
    end
    
end

%% save new job.mat (will have changed due to incidence change)
mises_directory = strrep(directory,'TURBOSTREAM','MISES');
cd(mises_directory)
% save new job
save([mises_directory 'job.mat'],'job');

%% change ises back to normal smoothing AND RERUN
display('******* Run at a reduced smoothing *******')
% Write ises file
fid = fopen([directory 'ises.mises'],'w');
fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');

% supersonic
fprintf(fid,'%s\n',' 15 4 3 17 6 !Global constraints');


fprintf(fid,'%6.5f %6.5f %6.5f %6.5f %6.5f ',[f.M f.Pin/f.Poin tand(f.Alpha) -0.66*c.m_chord f.Vtao1]);
fprintf(fid,'%s\n','|Minl p1/p01 Sinl Xinl v1/ao1');
fprintf(fid,'%6.5f %6.5f %6.5f %6.5f  ', [f.M2 f.Pout/f.Poin tand(f.Alpha2) 1.33*c.m_chord]);
fprintf(fid,'%s\n','|Mout p2/p01 Sout Xout');
fprintf(fid,'%s\n','0.00000     0.00000     1.39433     0.00000  | MFR  HWRATin GAMin DIFACTOR');
fprintf(fid,'%i ',round(f.Re));
fprintf(fid,'%s\n', '-4.00000 0.00000  0.26133 |Re -Turb');
if is_turb == 1
    fprintf(fid,'%s\n','0.001 0.001  |Xtr1 Xtr2');
else
    fprintf(fid,'%s\n','1.10 1.10  |Xtr1 Xtr2');
end
fprintf(fid,'%s\n','4     0.65000  -4.00000  0.00000   | ISMOM  MCRIT  MUCON  PLOSSIN');
fprintf(fid,'%s\n','0.00000  0.00000   | BVR1in  BVR2in');
fclose(fid);
% p=[];

if p.manual == 0;
    if run_single == 1
        
        [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_ss.py"','-echo');
        % run with multiple gridpar gain if it fails to converge
        [Polarx, Ises] = mis_read_polarx('mises',directory);
        if isstruct(Polarx) ~= 1
            n_type = 1;
            disp('*******Run again*********');
            while isstruct(Polarx) ~= 1 && n_type < n_type_max+1
                n_type = n_type + 1;
                disp(['****Create new gridpar ' num2str(n_type) '****']);
                mis_gridpar(directory,n_type);
                
                % subsonic
                [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_ss.py"','-echo');
                
                p.manual = 0;
                [Polarx, Ises] = mis_read_polarx('mises',directory);
            end
            if n_type> n_type_max
                disp('*******Have to run manually*********');
                while isstruct(Polarx) ~= 1
                    p.manual = 1;
                    Polarx = []; Polarx.p = {};
                    break
                end
            end
            n_type = 1; mis_gridpar(directory,n_type);
        end
        
        if exist('job','var') == 0
            [p,h] = mis_plot_section(directory,h,col,plot_stuff);
        else
            [p,h] = mis_plot_section(directory,h,col,plot_stuff,job);
        end
        if isempty(p) == 1
            p.manual = 1; % state it has not converged at reduced smoothing
            disp('not converged at reduced smoothing')
        else
            p.manual = 0;
        end
        
    else
        %% CREATE CODE TO GO THROUGH GRIDPAR AGAIN THERE IS SOME SMALL DIFFERENNCE IN THE SINGLE AND Ho-on FUNCTION  (in reading idat file?)
        [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMises_Ho-on_ss.py"','-echo');
        [Polarx, Ises] = mis_read_polarx('mises',directory);
        if isstruct(Polarx) ~= 1
            n_type = 1;
            disp('*******Run again*********');
            while isstruct(Polarx) ~= 1 && n_type < n_type_max+1
                n_type = n_type + 1;
                disp(['****Create new gridpar ' num2str(n_type) '****']);
                mis_gridpar(directory,n_type);
                
                % supersonic
                [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMises_Ho-on_ss.py"','-echo');
                
                p.manual = 0;
                [Polarx, Ises] = mis_read_polarx('mises',directory);
            end
            if n_type> n_type_max
                disp('*******Have to run manually*********');
                while isstruct(Polarx) ~= 1
                    p.manual = 1;
                    Polarx = []; Polarx.p = {};
                    break
                end
            end
            n_type = 1; mis_gridpar(directory,n_type);
        end
        
        % plot loss loop
        if plot_stuff == 1
            [Polar]=mis_read_polar('mises',directory);
            plot(Polar.ainl-abs(f.Alpha),Polar.yp,'Color',col);
            ylabel(['Y_{p}'],'Rotation',0)
            xlabel(['\Delta \alpha'])
        end
        
    end
    
    job.rjm_directory = job.rjm_directory_temp;
end

end


