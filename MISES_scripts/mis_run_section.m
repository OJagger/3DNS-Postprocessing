function [p,h,job] = mis_run_section(directory,h,col,plot_stuff,is_turb,run_single,job)
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
    % determine target Mach number triangles
    if isfield(job,'f') == 0
        [f]=mises_output_finput(job);
    else
        f=job.f;
        f.Re = 1*10^6;
        f.M = job.Mrel_inlet;
        f.AVDR = job.f.Ax;
    end
    if isfield(job,'psi_target') == 1
        psi_target = job.psi_target;
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
%     if job.Axt ~= 1
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
        
%     else
%         
%         % apply radial contraction from 25% chord upstream to 25% chord downstream
%         fprintf(fid,'%i\n',0);
%         fprintf(fid,'%6.5f %6.5f %6.5f\n',[min(mt(:,1))-0.5 c.r_le/c.tchord 1.00000]);
%         fprintf(fid,'%6.5f %6.5f %6.5f\n',[0.0-0.25*c.m_chord c.r_le/c.tchord 1.00000]);
%         fprintf(fid,'%6.5f %6.5f %6.5f\n',[0.0-0.25*c.m_chord c.r_le/c.tchord 1.00000]);
%         fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1))+0.25*c.m_chord c.r_te/c.tchord f.AVDR]);
%         fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1))+0.25*c.m_chord c.r_te/c.tchord f.AVDR]);
%         fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1))+0.5 c.r_te/c.tchord f.AVDR]);
%         fclose(fid);
%     end
end

%% Write ises file
fid = fopen(fullfile(directory, 'ises.mises'),'w');
fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
% subsonic
if f.M<1
    fprintf(fid,'%s\n',' 1 3 4 6 15 !Global constraints');
else
    % supersonic
    fprintf(fid,'%s\n',' 15 4 3 17 6 !Global constraints');
end

fprintf(fid,'%6.5f %6.5f %6.5f %6.5f %6.5f ',[f.M f.Pin/f.Poin tand(f.Alpha) -1.5*c.m_chord f.Vtao1]);
fprintf(fid,'%s\n','|Minl p1/p01 Sinl Xinl v1/ao1');
% fprintf(fid,'0.00000 0.00000 0.00000 %6.5f |Mout p2/p01 Sout Xout\n',1.33*m_chord);
fprintf(fid,'%6.5f %6.5f %6.5f %6.5f  ', [f.M2 f.Pout/f.Poin tand(f.Alpha2) 2.5*c.m_chord]);
fprintf(fid,'%s\n','|Mout p2/p01 Sout Xout');
fprintf(fid,'%s\n','0.00000     0.00000     1.39433     0.00000  | MFR  HWRATin GAMin DIFACTOR');
fprintf(fid,'%i ',round(f.Re));
fprintf(fid,'%s\n', '-4.00000 0.00000  0.26133 |Re -Turb');
if is_turb == 1
    fprintf(fid,'%s\n','0.001 0.001  |Xtr1 Xtr2');
else
    fprintf(fid,'%s\n','1.10 1.10  |Xtr1 Xtr2');
end
fprintf(fid,'%s\n','4     0.95000  -4.00000  0.00000   | ISMOM  MCRIT  MUCON  PLOSSIN');
fprintf(fid,'%s\n','0.00000  0.00000   | BVR1in  BVR2in');
fprintf(fid,'%s\n','5.6  1.0  1.0  | SCC SCP SCD');
fprintf(fid,'%s\n','-0.071652  1.0  | XSHOCK FCTSHK');
fclose(fid);


%% create mesh
n_type = 1;
n_type_max = 7;
mis_gridpar(directory,n_type);

% Change to MISES directory and set environment variables
cd(directory)
setenv('GFORTRAN_STDIN_UNIT', '5')
setenv('GFORTRAN_STDOUT_UNIT', '6')
setenv('GFORTRAN_STDERR_UNIT', '0')
setenv('LD_LIBRARY_PATH','/opt/intel/fce/10.1.015/lib:/home/dl467/lib:/home/dl467/lib')

%% Run MISES solver

% Run MISES and change directory back
% subsonic
[~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle.py"','-echo');

cd(directory)
p.manual = 0;

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
                
                % subsonic
                [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle.py"','-echo');
                
                p.manual = 0;
                [Polarx, Ises] = mis_read_polarx('mises',directory);
            else
                disp('*******Have to run manually*********');
                while isstruct(Polarx) ~= 1
                    p.manual = 1;
                    Polarx = [];
                    Polarx.p = {};
                    break
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
            % subsonic
            [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle.py"','-echo');
            p.manual = 0;
            [Polarx, Ises] = mis_read_polarx('mises',directory);
        else
            disp('*******Have to run manually*********');
            while isstruct(Polarx) ~= 1
                p.manual = 1;
                Polarx = [];
                Polarx.p = {};
                break
            end
        end
    end
    n_type = 1;
    % 	return
end

%% Read in the grid coodinates
if exist([directory 'idat.mises_01'],'file') ~= 0
    Idat = mis_read_idat('mises_01',directory);
else
    Idat = mis_read_idat('mises',directory);
end


%% correct for incidence
if exist('psi_target','var')== 1
    % Calculate the relative angle from the leading edge centre
    mt_stag = [Idat.x(Idat.ninl(2),end) Idat.y(Idat.ninl(2),end)];
    phi_stag = atan2(c.mt_le_cen(2) - mt_stag(2),c.mt_le_cen(1) - mt_stag(1)) * 360 / (2*pi);
    mt = c.mt(c.mt(:,1) < c.mt_le_cen(1),:);
    phi = atan2(c.mt_le_cen(2) - mt(:,2),c.mt_le_cen(1) - mt(:,1)) * 360 / (2*pi);
    
    % Calculate the angle the stagnation point makes with the leading edge
    [phi,i] = unique(phi); mt = mt(i,:);
    dtdx = grad_mg(mt(:,1),mt(:,2));
    psi = atand(dtdx);
    
    % Correct the angle by the leading edge metal angle
    p.psi_stag = 90 - (c.chi_le - interp1(phi,psi,phi_stag,'pchip'));
    psi_stag = p.psi_stag;
    
    if isnan(psi_stag) == 1
        psi_stag = psi_target + 10;
        disp('***** psi solution is nan ******')
    end
    
    % iteration number
    ni = 1;
    n_psi_max = 5;
    n_type_max = 7;
    beta=0.05;
    p.manual = 0;
    while abs(psi_stag-psi_target)>5.0 && ni < n_psi_max 
%         && isstruct(Polarx) == 1 && p.manual == 0
        
        disp(['***** Difference Dpsi = ' num2str(psi_stag-psi_target) ' *******'])
        % run design and correct for incidence
        [~,~]=system(['rm ' directory 'polarx.mises'],'-echo');
        if ni==1
            
            dchi_le = (psi_stag-psi_target)*beta
            if abs(dchi_le) > 1
                dchi_le = 1*sign(dchi_le);
            end
            chi_le = job.chi_le - dchi_le;
            job.chi_le_0 = job.chi_le;
            % 			job.psi_0 = delta.psi;
            job.psi_0 = psi_stag;
            job.chi_le=chi_le;
            
            % keep total camber constant
            job.chi_te_pre = job.chi_te;
            dchi_te = - dchi_le;
            job.chi_te = job.chi_te + dchi_te;
            display(['Total camber is :' num2str(c.chi_le-c.chi_te)])
            [mt,c] = mises_create_mt(job,1);
            % Write section file
            save([directory 'section.mat'],'c','f');
            % subsonic
            [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle.py"','-echo');
            
            p.manual = 0;
        else
            
            if abs(dchi_le) > 1
                dchi_le = 1*sign(dchi_le);
            end
            chi_le = job.chi_le - dchi_le;
            job.chi_le_0 = job.chi_le;
            % 			job.psi_0 = delta.psi;
            job.psi_0 = psi_stag;
            job.chi_le = chi_le;
            
            %keep total camber constant
            job.chi_te_pre = job.chi_te;
            dchi_te = - dchi_le;
            job.chi_te = job.chi_te + dchi_te;
            display(['Total camber is :' num2str(c.chi_le-c.chi_te)])
            
            [mt,c] = mises_create_mt(job,1);
            
            % store paramaterised blade parameters and Mach number triangles in section.mat
            save([directory 'section.mat'],'c','f');
            % subsonic
            [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle.py"','-echo');
            
            p.manual = 0;
            
        end
        
        ni = ni+1;
        
        % Read in flow file and re-run if not converged with a new mesh
        if exist([directory 'polarx.mises'],'file') ~= 0
            [Polarx, Ises] = mis_read_polarx('mises',directory);
            if isstruct(Polarx) ~= 1
                disp('*******Run again*********');
                while isstruct(Polarx) ~= 1 && n_type < n_type_max + 1
                    n_type = n_type + 1;
                    disp(['****Create new gridpar ' num2str(n_type) '****']);
                    mis_gridpar(directory,n_type);
                    
                    % subsonic
                    [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle.py"','-echo');
                    
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
        else
            disp('File Not Found')
            disp('*******Run again*********');
            Polarx = [];
            while isstruct(Polarx) ~= 1 && n_type < n_type_max + 1
                n_type = n_type + 1;
                disp(['****Create new gridpar ' num2str(n_type) '****']);
                mis_gridpar(directory,n_type);
                
                % subsonic
                [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle.py"','-echo');
                
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
        
        % Read in the grid coodinates
        if exist([directory 'idat.mises_01'],'file') ~= 0
            Idat = mis_read_idat('mises_01',directory);
        else
            Idat = mis_read_idat('mises',directory);
        end
        
        % Calculate the relative angle from the leading edge centre
        mt_stag = [Idat.x(Idat.ninl(2),end) Idat.y(Idat.ninl(2),end)];
        phi_stag = atan2(c.mt_le_cen(2) - mt_stag(2),c.mt_le_cen(1) - mt_stag(1)) * 360 / (2*pi);
        mt = c.mt(c.mt(:,1) < c.mt_le_cen(1),:);
        phi = atan2(c.mt_le_cen(2) - mt(:,2),c.mt_le_cen(1) - mt(:,1)) * 360 / (2*pi);
        
        % Calculate the angle the stagnation point makes with the leading edge
        [phi,i] = unique(phi); mt = mt(i,:);
        dtdx = grad_mg(mt(:,1),mt(:,2));
        psi = atand(dtdx);
        
        % Correct the angle by the leading edge metal angle
        p.psi_stag = 90 - (c.chi_le - interp1(phi,psi,phi_stag,'pchip'));
        psi_stag = p.psi_stag;
        if isnan(psi_stag) == 1
            psi_stag = psi_target + 0.5* dchi_le ./ beta;
            dchi_le = (psi_stag-psi_target)*beta
            job.chi_le = job.chi_le_0;
            disp('***** psi solution is nan - half the change ******')
        else
            %provide a limiter to increase in kappa_in
            dchi_le = (job.chi_le-job.chi_le_0)/(psi_stag-job.psi_0)*(psi_stag-psi_target)
        end
    end
end

%% save new job.mat (will have changed due to incidence change)
mises_directory = strrep(directory,'TURBOSTREAM','MISES');
cd(mises_directory)
% save new job
save([mises_directory 'job.mat'],'job');

%% change ises back to normal smoothing and run
display('******* Run at a reduced smoothing *******')
% Write ises file
fid = fopen([directory 'ises.mises'],'w');
fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');

% subsonic
fprintf(fid,'%s\n',' 1 3 4 6 15 !Global constraints');



fprintf(fid,'%6.5f %6.5f %6.5f %6.5f %6.5f ',[f.M f.Pin/f.Poin tand(f.Alpha) -0.66*c.m_chord f.Vtao1]);
fprintf(fid,'%s\n','|Minl p1/p01 Sinl Xinl v1/ao1');
fprintf(fid,'%6.5f %6.5f %6.5f %6.5f  ', [f.M2 f.Pout/f.Poin tand(f.Alpha2) 1.33*c.m_chord]);
fprintf(fid,'%s\n','|Mout p2/p01 Sout Xout');
fprintf(fid,'%s\n','0.00000     0.00000     1.39433     0.00000  | MFR  HWRATin GAMin DIFACTOR');
fprintf(fid,'%i ',round(f.Re));
fprintf(fid,'%s\n', '-4.00000 0.00000  0.26133 |Re -Turb');
if is_turb == 1
    fprintf(fid,'%s\n','0.001 0.01  |Xtr1 Xtr2');
else
    fprintf(fid,'%s\n','1.10 1.10  |Xtr1 Xtr2');
end
fprintf(fid,'%s\n','4     0.95000  -1.00000  0.00000   | ISMOM  MCRIT  MUCON  PLOSSIN');
fprintf(fid,'%s\n','0.00000  0.00000   | BVR1in  BVR2in');
fclose(fid);


if p.manual == 0;
    if run_single == 1
        
        % subsonic
        [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle.py"','-echo');
        % run gridpar again if it fails to  converge
        [Polarx, Ises] = mis_read_polarx('mises',directory);
        if isstruct(Polarx) ~= 1
            n_type = 1;
            disp('*******Run again*********');
            while isstruct(Polarx) ~= 1 && n_type < n_type_max+1
                n_type = n_type + 1;
                disp(['****Create new gridpar ' num2str(n_type) '****']);
                mis_gridpar(directory,n_type);
                
                % subsonic
                [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMises_Ho-on.py"','-echo');
                
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
        %% CREATE CODE TO GO THROUGH GRIDPAR AGAIN THERE IS SOME SMALL DIFFERENNCE IN THE SINGLE AND Ho-on FUNCTION (in reading idat file?)
        [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMises_Ho-on.py"','-echo');
        [Polarx, Ises] = mis_read_polarx('mises',directory);
        if isstruct(Polarx) ~= 1
            n_type = 1;
            disp('*******Run again*********');
            while isstruct(Polarx) ~= 1 && n_type < n_type_max+1
                n_type = n_type + 1;
                disp(['****Create new gridpar ' num2str(n_type) '****']);
                mis_gridpar(directory,n_type);
                
                % subsonic
                [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMises_Ho-on.py"','-echo');
                
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
        
        if plot_stuff == 1
            [Polar]=mis_read_polar('mises',directory);
            plot(Polar.ainl-abs(f.Alpha),Polar.yp,'Color',col);
            ylabel(['Y_{p}'],'Rotation',0)
            xlabel(['\Delta \alpha'])
        end
        
    end
end

job.rjm_directory = job.rjm_directory_temp;

end

