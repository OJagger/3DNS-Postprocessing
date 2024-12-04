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

% Prepare defaults
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

% 
% exist('job','var')


if exist('job','var') == 1
% 	[c]=mises_output_xrrt(job);
	if isfield(job,'f') == 0
		[f]=mises_output_finput(job);
	else
		f=job.f;
		f.Re = 1*10^6;
		f.M = job.Mrel_inlet;
		f.AVDR = job.f.Ax;
    end
end


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

save([directory 'job.mat'],'job');


if f.Alpha < 0; f.Alpha = -f.Alpha; end;
if f.Alpha2 < 0; f.Alpha2 = -f.Alpha2; end;

%% Prepare coordinates for MISES

[mt,c] = mises_create_mt(job,1);
% Maximum MISES length
m_chord = max(mt(:,1));

% determine throat location

%% Write input files
% Write section file
save([directory 'section.mat'],'c','f');

%%


% Write blade file
fid = fopen([directory 'blade.mises'],'w');
fprintf(fid,'%s\n','MISES_Section');
fprintf(fid,'%6.5f %6.5f %3.2f %3.2f %7.6f\n',...
    [tand(c.chi_le) tand(c.chi_te) m_chord m_chord 2*pi / c.N]);
% mt(:,2) = - mt(:,2);
for i = 1:size(mt,1)
    fprintf(fid,'%10.12f %10.12f\n',mt(i,:));
end
fclose(fid);

%% Write stream file
fid = fopen([directory 'stream.mises'],'w');

% apply radial contraction from 25% chord upstream to 25% chord downstream
fprintf(fid,'%i\n',0);
fprintf(fid,'%6.5f %6.5f %6.5f\n',[min(mt(:,1))-0.5 c.r_le/c.tchord 1.00000]);
fprintf(fid,'%6.5f %6.5f %6.5f\n',[0.0-0.25*c.m_chord c.r_le/c.tchord 1.00000]);
fprintf(fid,'%6.5f %6.5f %6.5f\n',[0.0-0.25*c.m_chord  c.r_le/c.tchord 1.00000]);
fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1))+0.25*c.m_chord c.r_te/c.tchord f.AVDR]);
fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1))+0.25*c.m_chord c.r_te/c.tchord f.AVDR]);
fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1))+0.5 c.r_te/c.tchord f.AVDR]);
%% defines the inlet (do not change unless radial contraction starts earlier than LE)
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[min(mt(:,1))-0.5 c.r_le/c.tchord 1.00000]);
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[min(mt(:,1))-0.5 c.r_le/c.tchord 1.00000]);
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[0.0 c.r_le/c.tchord 1.00000]);
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[0.0 c.r_le/c.tchord 1.00000]);
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[0.0 c.r_le/c.tchord 1.00000]);
%% changed to start radial contraction from 50% DONT FORGET TO CHANGE
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[-0.20*max(mt(:,1)) c.r_le/c.tchord 1.00000]);
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[-0.20*max(mt(:,1)) c.r_le/c.tchord 1.00000]);
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[-0.0*max(mt(:,1)) c.r_le/c.tchord 1.00000]);
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[0.80*max(mt(:,1)) c.r_le/c.tchord 1-(1-f.AVDR)*0.25]);
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1))+0.40*max(mt(:,1)) c.r_te/c.tchord f.AVDR]);
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1))+0.40*max(mt(:,1)) c.r_te/c.tchord f.AVDR]);
%% changed radial contraction up to throat and back again
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[c.xs_sb*max(mt(:,1)) c.r_le/c.tchord f.AVDR]);
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[c.xs_sb*max(mt(:,1)) c.r_le/c.tchord f.AVDR]);
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[c.xs_sb*max(mt(:,1)) c.r_le/c.tchord f.AVDR]);
% % fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1))+0.40*max(mt(:,1)) c.r_te/c.tchord 1]);
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1)) c.r_te/c.tchord 1]);
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1)) c.r_te/c.tchord 1]);
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1))+0.50 c.r_te/c.tchord 1]);
fclose(fid);

% Write ises file
fid = fopen([directory 'ises.mises'],'w');
% fprintf(fid,'%s\n','1 2 5  !Global variables');
% fprintf(fid,'%s\n','1 3 4  !Global constraints');
% fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
% fprintf(fid,'%s\n',' 3 4 6 9 17 !Global constraints');
fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
%% subsonic
% fprintf(fid,'%s\n',' 1 3 4 6 15 !Global constraints');
%% supersonic
fprintf(fid,'%s\n',' 15 4 3 17 6 !Global constraints');

fprintf(fid,'%6.5f %6.5f %6.5f %6.5f %6.5f ',[f.M f.Pin/f.Poin tand(f.Alpha) -0.66*m_chord f.Vtao1]);
fprintf(fid,'%s\n','|Minl p1/p01 Sinl Xinl v1/ao1');
% fprintf(fid,'0.00000 0.00000 0.00000 %6.5f |Mout p2/p01 Sout Xout\n',1.33*m_chord);
fprintf(fid,'%6.5f %6.5f %6.5f %6.5f  ', [f.M2 f.Pout/f.Poin tand(f.Alpha2) 1.33*m_chord]);
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


%% Run MISES solver and plot results

% Loop over multiple grid sizes until converged
% for n = 78:-1:70

n=100;

% % Write gridpar file
% fid = fopen([directory 'gridpar.mises'],'w');
% fprintf(fid,'%s\n','T T');
% fprintf(fid,'%s\n','60   45');
% fprintf(fid,'%s\n','24');
% fprintf(fid,'%s\n','1.500000');
% fprintf(fid,'%s\n','0.800000');
% fprintf(fid,'%s\n',[num2str(n) '    0.10000    0.900000    1.000000']);
% fprintf(fid,'%s\n','1.000000    1.000000    0.000000    1.000000    1.000000    0.000000');
% fclose(fid);
n_type_max = 7;
n_type = 1; % n=100
mis_gridpar(directory,n_type);

% Change to MISES directory and set environment variables
cd(directory)
setenv('GFORTRAN_STDIN_UNIT', '5')
setenv('GFORTRAN_STDOUT_UNIT', '6')
setenv('GFORTRAN_STDERR_UNIT', '0')
setenv('LD_LIBRARY_PATH','/opt/intel/fce/10.1.015/lib:/home/dl467/lib:/home/dl467/lib')

% Run MISES and change directory back
if f.M<1
%% subsonic
[~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle.py"','-echo');
else
%% supersonic
[~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_ss.py"','-echo');
end
p.manual = 0;
cd(directory)

% Read in flow file and re-run if not converged
if exist([directory 'polarx.mises'],'file') ~= 0
	[Polarx, Ises] = mis_read_polarx('mises',directory);
	if isstruct(Polarx) ~= 1
		disp('*******Run again*********');		
		while isstruct(Polarx) ~= 1 
			n_type = n_type + 1
			
			if n_type<=n_type_max
			disp(['****Create new gridpar ' num2str(n_type) '****']);
			mis_gridpar(directory,n_type);
			if f.M<1
			%% subsonic
			[~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle.py"','-echo');
			else
			%% supersonic
			[~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_ss.py"','-echo');
			end
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
						[~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_manual_ss.py"','-echo');
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
		n_type = n_type + 1
		if n_type<=n_type_max+1
		disp(['****Create new gridpar ' num2str(n_type) '****']);
		mis_gridpar(directory,n_type);

		if f.M<1
		%% subsonic
		[~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle.py"','-echo');
		else
		%% supersonic
		[~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_ss.py"','-echo');
		end
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
					[~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_manual_ss.py"','-echo');
				end
				[Polarx, Ises] = mis_read_polarx('mises',directory);
			end
		end
	end
	n_type = 1;
% 	p = [];
% 	return
end

% Read in the grid coodinates
if exist([directory 'idat.mises_01'],'file') ~= 0
    Idat = mis_read_idat('mises_01',directory);
else
    Idat = mis_read_idat('mises',directory);
end


%% correct for incidence
% if exist('psi_target','var')== 1
% 	
% 	% Calculate the angle the stagnation point makes with the leading edge
% 	[phi,i] = unique(phi); mt = mt(i,:);
% 	dtdx = grad_mg(mt(:,1),mt(:,2));
% 	psi = atand(dtdx);
% 	
% 	% Correct the angle by the leading edge metal angle
% 	p.psi_stag = 90 - (c.chi_le - interp1(phi,psi,phi_stag,'pchip'));
% 	psi_stag = p.psi_stag;

	[Polarx, Ises] = mis_read_polarx('mises',directory);
	binl_target = Ises.binl;
	if p.manual == 1
		disp('*****Cannot run initial iteration*****')
		binl = binl_target;
	else
	binl = Polarx.binl;
    % Calculate the relative angle from the leading edge centre
	mt_stag = [Idat.x(Idat.ninl(2),end) Idat.y(Idat.ninl(2),end)];
	phi_stag = atan2(c.mt_le_cen(2) - mt_stag(2),c.mt_le_cen(1) - mt_stag(1)) * 360 / (2*pi);
	mt = c.mt(c.mt(:,1) < c.mt_le_cen(1),:);
	phi = atan2(c.mt_le_cen(2) - mt(:,2),c.mt_le_cen(1) - mt(:,1)) * 360 / (2*pi);
	end
	
	% iteration number
	ni = 1;
	n_psi_max = 10;
	n_type_max = 7;
	beta=0.05;
	while abs(binl_target-binl)>0.15 && ni < n_psi_max && isstruct(Polarx) == 1 && p.manual == 0

		[~,~]=system(['rm ' directory 'polarx.mises'],'-echo');
% 		if ni==1
			
			dinlet_cam = binl_target-binl;
% 			psi_stag
% 			psi_target
% 			dinlet_cam = (psi_stag-psi_target)*beta
			if abs(dinlet_cam) > 1
				dinlet_cam = 1*sign(dinlet_cam);
			end
			inlet_cam = job.inlet_cam + dinlet_cam;
			job.inlet_cam_0 = job.inlet_cam;
			binl_0 = binl;
% % 			job.psi_0 = delta.psi;
% 			job.psi_0 = psi_stag;
			job.inlet_cam=inlet_cam;
% 
% 			% change total camber to keep same exit angle
% 			job.tot_cam_pre = job.tot_cam;
% 			job.tot_cam = job.tot_cam + dinlet_cam;
			% change percentage ss camber to keep ss section the
			% same
% 			job.percentage_cam = double(double(job.percentage_cam) * double(job.tot_cam_pre)/double(job.tot_cam));
			[mt,c] = mises_create_mt(job,1);
			% Write section file
			save([directory 'section.mat'],'c','f');
			%% subsonic
			if f.M<1
			[~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle.py"','-echo');
			else
			%% supersonic
			[~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_ss.py"','-echo');
			end
% 		else
% 			%provide a limiter to increase in kappa_in
% %             psi_target
% % 			psi_stag
% % 			job.psi_0
% % 			dinlet_cam = (job.inlet_cam-job.inlet_cam_0)/(delta.psi-job.psi_0)*(delta.psi-job.psi_target)
% 
%             dinlet_cam = (job.inlet_cam-job.inlet_cam_0)/(binl-binl_0)*(binl-binl_target)
% 			if abs(dinlet_cam) > 1
% 				dinlet_cam = 1*sign(dinlet_cam);
% 			end
% 			inlet_cam = job.inlet_cam + dinlet_cam;
% 			job.inlet_cam_0 = job.inlet_cam;
% 			binl_0 = binl;
% % 			job.psi_0 = delta.psi;
% %             job.psi_0 = psi_stag;
% % 			job.inlet_cam = inlet_cam;
% % 			job.tot_cam_pre = job.tot_cam;
% % 			% 		% change total camber to keep same exit angle
% 			job.tot_cam = job.tot_cam + dinlet_cam;
% 			% change percentage ss camber to keep ss section the
% 			% same
% % 			job.percentage_cam = double(double(job.percentage_cam) * double(job.tot_cam_pre)/double(job.tot_cam));
% 			[mt,c] = mises_create_mt(job,1);
% 			% Write section file
% 			save([directory 'section.mat'],'c','f');
% 			if f.M<1
% 			%% subsonic
% 			[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle.py"','-echo');
% 			else
% 			%% supersonic
% 			[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle_ss.py"','-echo');
% 			end
% 
% 
% 		end
		
		ni = ni+1;

			% Read in flow file and re-run if not converged
			if exist([directory 'polarx.mises'],'file') ~= 0
				[Polarx, Ises] = mis_read_polarx('mises',directory);
				if isstruct(Polarx) ~= 1
					disp('*******Run again*********');
					while isstruct(Polarx) ~= 1 && n_type < n_type_max+1
						n_type = n_type + 1;
						disp(['****Create new gridpar ' num2str(n_type) '****']);
						mis_gridpar(directory,n_type);
						if f.M<1
						%% subsonic
						[~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle.py"','-echo');
						else
						%% supersonic
						[~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_ss.py"','-echo');
						[Polarx, Ises] = mis_read_polarx('mises',directory);
						end
						p.manual = 0;
					end
					if n_type> n_type_max
						disp('*******Have to run manually*********');
						while isstruct(Polarx) ~= 1
							if f.M<1
								[~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_manual.py"','-echo');
							else
								p.manual = 1;
								Polarx = []; Polarx.p = {};
								break
								[~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_manual_ss.py"','-echo');
							end
							
							[Polarx, Ises] = mis_read_polarx('mises',directory);
						end
					end
					n_type = 1; mis_gridpar(directory,n_type);
				end
			else
				disp('File Not Found')
				disp('*******Run again*********');
				Polarx = [];
				while isstruct(Polarx) ~= 1 && n_type < n_type_max+1
					n_type = n_type + 1;
					disp(['****Create new gridpar ' num2str(n_type) '****']);
					mis_gridpar(directory,n_type);
					if f.M<1
					%% subsonic
					[~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle.py"','-echo');
					else
					%% supersonic
					[~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_ss.py"','-echo');
					[Polarx, Ises] = mis_read_polarx('mises',directory);
					end
				end
				p.manual = 0;
				if n_type> n_type_max
					disp('*******Have to run manually*********');
					while isstruct(Polarx) ~= 1
						if f.M<1
							[~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_manual.py"','-echo');
						else
						p.manual = 1;
						Polarx = []; Polarx.p = {};
						break
						[~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMisesSingle_manual_ss.py"','-echo');
						end
						[Polarx, Ises] = mis_read_polarx('mises',directory);
					end
				end
				n_type = 1; mis_gridpar(directory,n_type);
% 				p = [];
				% 			return
			end
		
		% Read in the grid coodinates
		if exist([directory 'idat.mises_01'],'file') ~= 0
			Idat = mis_read_idat('mises_01',directory);
		else
			Idat = mis_read_idat('mises',directory);
		end
		
		% 		% Calculate the relative angle from the leading edge centre
		% 		mt_stag = [Idat.x(Idat.ninl(2),end) Idat.y(Idat.ninl(2),end)];
		% 		phi_stag = atan2(c.mt_le_cen(2) - mt_stag(2),c.mt_le_cen(1) - mt_stag(1)) * 360 / (2*pi);
		% 		mt = c.mt(c.mt(:,1) < c.mt_le_cen(1),:);
		% 		phi = atan2(c.mt_le_cen(2) - mt(:,2),c.mt_le_cen(1) - mt(:,1)) * 360 / (2*pi);
		%
		% 		% Calculate the angle the stagnation point makes with the leading edge
		% 		[phi,i] = unique(phi); mt = mt(i,:);
		% 		dtdx = grad_mg(mt(:,1),mt(:,2));
		% 		psi = atand(dtdx);
		%
		% 		% Correct the angle by the leading edge metal angle
		% 		p.psi_stag = 90 - (c.chi_le - interp1(phi,psi,phi_stag,'pchip'));
		% 		psi_stag = p.psi_stag;
		if exist([directory 'polarx.mises'],'file') ~= 0 && p.manual == 0
			[Polarx, Ises] = mis_read_polarx('mises',directory);
			binl_target = Ises.binl
			binl = Polarx.binl
		end
% 	end
end

%% save new job.mat (will have changed due to incidence change)
mises_directory = strrep(directory,'TURBOSTREAM','MISES');
cd(mises_directory)
% save new job
save([mises_directory 'job.mat'],'job');

%% change ises back to normal smoothing
% Write ises file
fid = fopen([directory 'ises.mises'],'w');
% fprintf(fid,'%s\n','1 2 5  !Global variables');
% fprintf(fid,'%s\n','1 3 4  !Global constraints');
% fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
% fprintf(fid,'%s\n',' 3 4 6 9 17 !Global constraints');
fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
if f.M<1
%% subsonic
fprintf(fid,'%s\n',' 1 3 4 6 15 !Global constraints');
else
%% supersonic
fprintf(fid,'%s\n',' 15 4 3 17 6 !Global constraints');
end

fprintf(fid,'%6.5f %6.5f %6.5f %6.5f %6.5f ',[f.M f.Pin/f.Poin tand(f.Alpha) -0.66*m_chord f.Vtao1]);
fprintf(fid,'%s\n','|Minl p1/p01 Sinl Xinl v1/ao1');
% fprintf(fid,'0.00000 0.00000 0.00000 %6.5f |Mout p2/p01 Sout Xout\n',1.33*m_chord);
fprintf(fid,'%6.5f %6.5f %6.5f %6.5f  ', [f.M2 f.Pout/f.Poin tand(f.Alpha2) 1.33*m_chord]);
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
        % 	if f.M<1
        % 		%% subsonic
        % 		[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle.py"','-echo');
        % 	else
        % 		%% supersonic
        % 		[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle_ss.py"','-echo');
        % 		[Polarx, Ises] = mis_read_polarx('mises',directory);
        % 	end
        if exist('job','var') == 0
            [p,h] = mis_plot_section(directory,h,col,plot_stuff);
        else
            [p,h] = mis_plot_section(directory,h,col,plot_stuff,job);
        end
        p.manual = 0;
    else
        if f.M<1
            [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMises_Ho-on.py"','-echo');
        else
            [~,~] = system('bash --login -c "python /mnt/Disk1/dl467/bin/CalculateBladeMises_Ho-on_ss.py"','-echo');
        end
        if plot_stuff == 1
            [Polar]=mis_read_polar('mises',directory);
            plot(Polar.ainl-abs(f.Alpha),Polar.yp,'Color',col);
            ylabel(['Y_{p}'],'Rotation',0)
            xlabel(['\Delta \alpha'])
            p=[];
        end
        p.manual = 0;
    end

job.rjm_directory = job.rjm_directory_temp;
end

end


