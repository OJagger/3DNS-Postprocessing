function [p,h] = mis_run_section_manufacture(directory,h,col,plot_stuff,is_turb,run_single,job,process)
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
% process can either be string: 'forge' or 'ecm'

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

if ischar(process)==1
dr.mf  = '/home/dl467/Documents/Manufacturing';
cd(dr.mf);
load([dr.mf '/' 'dev_' process '2.mat']);
% load([dr.mf '/' 'dev_' process '3.mat']);
end

if exist('prior_solution','var') == 0
	prior_solution = 1;
end

if exist('job','var') == 1
	[c]=mises_output_xrrt(job);
	[f]=mises_output_finput(job);
	g=ts_read_hdf5([job.rjm_directory job.outname]);
	[~,~,delta] = ts_incidence_2D(g,0);
	psi_target = delta.psi
end

directory = strrep(directory,'TURBOSTREAM','MISES');

% Make directory if required
if exist(directory,'dir') == 0
    mkdir(directory);
end

% if prior_solution == 0
% 
% % Delete old output files from directory
% % if isempty(dir([directory 'polar*'])) == 0; delete([directory 'polar*']); end;
% % if isempty(dir([directory 'output*'])) == 0; delete([directory 'output*']); end;
% % if isempty(dir([directory 'idat*'])) == 0; delete([directory 'idat*']); end;
% % if isempty(dir([directory 'spec*'])) == 0; delete([directory 'spec*']); end;
% 
% %% Prepare coordinates for MISES
% if f.Alpha < 0; f.Alpha = -f.Alpha; end;
% if f.Alpha2 < 0; f.Alpha2 = -f.Alpha2; end;
% 
% % if exist job.var use that to define xrt coordinates
% if exist('job','var')==1
% 	write_file = 1;
% 	[mt,c] = mises_create_mt(job,write_file);
% 	% Write section file
%     save([directory 'section.mat'],'c','f');
% end
% 
% %%
% 
% %% Write input files
% 
% % Write section file
% save([directory 'section.mat'],'c','f');
% 
% % Maximum MISES length
% m_chord = max(mt(:,1));
% r_le = c.r_le;
% r_te = c.r_te;
% 
% % Write blade file
% fid = fopen([directory 'blade.mises'],'w');
% fprintf(fid,'%s\n','MISES_Section');
% fprintf(fid,'%6.5f %6.5f %3.2f %3.2f %7.6f\n',...
%     [tand(c.chi_le) tand(c.chi_te) m_chord m_chord 2*pi / c.N]);
% % mt(:,2) = - mt(:,2);
% for i = 1:size(mt,1)
%     fprintf(fid,'%10.12f %10.12f\n',mt(i,:));
% end
% fclose(fid);
% 
% % Write stream file
% fid = fopen([directory 'stream.mises'],'w');
% fprintf(fid,'%i\n',0);
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[min(mt(:,1))-0.5 r_le/c.tchord 1.00000]);
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[0.0 r_le/c.tchord 1.00000]);
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[0.0 r_le/c.tchord 1.00000]);
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1)) r_te/c.tchord f.AVDR]);
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1)) r_te/c.tchord f.AVDR]);
% fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1))+0.5 r_te/c.tchord f.AVDR]);
% fclose(fid);
% 
% % Write ises file
% fid = fopen([directory 'ises.mises'],'w');
% % fprintf(fid,'%s\n','1 2 5  !Global variables');
% % fprintf(fid,'%s\n','1 3 4  !Global constraints');
% % fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
% % fprintf(fid,'%s\n',' 3 4 6 9 17 !Global constraints');
% fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
% fprintf(fid,'%s\n',' 1 3 4 6 15 !Global constraints');
% 
% fprintf(fid,'%6.5f %6.5f %6.5f %6.5f %6.5f ',[f.M f.Pin/f.Poin tand(f.Alpha) -0.66*m_chord f.Vtao1]);
% fprintf(fid,'%s\n','|Minl p1/p01 Sinl Xinl v1/ao1');
% % fprintf(fid,'0.00000 0.00000 0.00000 %6.5f |Mout p2/p01 Sout Xout\n',1.33*m_chord);
% fprintf(fid,'%6.5f %6.5f %6.5f %6.5f  ', [f.M2 f.Pout/f.Poin tand(f.Alpha2) 1.33*m_chord]);
% fprintf(fid,'%s\n','|Mout p2/p01 Sout Xout');
% fprintf(fid,'%s\n','0.00000     0.00000     1.39433     0.00000  | MFR  HWRATin GAMin DIFACTOR');
% fprintf(fid,'%i ',round(f.Re));
% fprintf(fid,'%s\n', '-4.00000 0.00000  0.26133 |Re -Turb');
% if is_turb == 1
%     fprintf(fid,'%s\n','0.001 0.01  |Xtr1 Xtr2');
% else
%     fprintf(fid,'%s\n','1.10 1.10  |Xtr1 Xtr2');
% end
% fprintf(fid,'%s\n','4     0.95000  -4.00000  0.00000   | ISMOM  MCRIT  MUCON  PLOSSIN');
% fprintf(fid,'%s\n','0.00000  0.00000   | BVR1in  BVR2in');
% fclose(fid);
% 
% 
% %% Run MISES solver and plot results
% 
% % Loop over multiple grid sizes until converged
% % for n = 78:-1:70
% 
% n=100;
% 
% % Write gridpar file
% % fid = fopen([directory 'gridpar.mises'],'w');
% % fprintf(fid,'%s\n','T T');
% % fprintf(fid,'%s\n','60   45');
% % fprintf(fid,'%s\n','24');
% % fprintf(fid,'%s\n','1.500000');
% % fprintf(fid,'%s\n','0.800000');
% % fprintf(fid,'%s\n',[num2str(n) '    0.10000    0.900000    1.000000']);
% % fprintf(fid,'%s\n','1.000000    1.000000    0.000000    1.000000    1.000000    0.000000');
% % fclose(fid);
% n_type = 1; % n=100
% mis_gridpar(directory,n_type);
% 
% % Change to MISES directory and set environment variables
% cd(directory)
% setenv('GFORTRAN_STDIN_UNIT', '5')
% setenv('GFORTRAN_STDOUT_UNIT', '6')
% setenv('GFORTRAN_STDERR_UNIT', '0')
% setenv('LD_LIBRARY_PATH','/opt/intel/fce/10.1.015/lib:/home/dl467/lib:/home/dl467/lib')
% 
% % Run MISES and change directory back
% [~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle.py"','-echo');
% cd(directory)
% 
% % Read in flow file and re-run if not converged
% if exist([directory 'polarx.mises'],'file') ~= 0
% 	[Polarx, Ises] = mis_read_polarx('mises',directory);
% 	if isstruct(Polarx) ~= 1
% 		disp('*******Run again*********');
% 		while isstruct(Polarx) ~= 1 && n_type < 4
% 			n_type = n_type + 1;
% 			disp(['****Create new gridpar ' num2str(n_type) '****']);
% 			mis_gridpar(directory,n_type);
% 			[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle.py"','-echo');
% 			[Polarx, Ises] = mis_read_polarx('mises',directory);
% 		end
% 		n_type = 1; mis_gridpar(directory,n_type);
% 	end
% else
% 	disp('File Not Found')
% 	disp('*******Run again*********');
% 	Polarx = [];
% 	while isstruct(Polarx) ~= 1 && n_type < 4
% 		n_type = n_type + 1;
% 		disp(['****Create new gridpar ' num2str(n_type) '****']);
% 		mis_gridpar(directory,n_type);
% 		[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle.py"','-echo');
% 		[Polarx, Ises] = mis_read_polarx('mises',directory);
% 	end
% 	n_type = 1; mis_gridpar(directory,n_type);
% 	p = [];
% 	% 			return
% end
% 
% % Read in the grid coodinates
% if exist([directory 'idat.mises_01'],'file') ~= 0
%     Idat = mis_read_idat('mises_01',directory);
% else
%     Idat = mis_read_idat('mises',directory);
% end
% 
% 
% %% correct for incidence
% if exist('psi_target','var')== 1
% 	% Calculate the relative angle from the leading edge centre
% 	mt_stag = [Idat.x(Idat.ninl(2),end) Idat.y(Idat.ninl(2),end)];
% 	phi_stag = atan2(c.mt_le_cen(2) - mt_stag(2),c.mt_le_cen(1) - mt_stag(1)) * 360 / (2*pi);
% 	mt = c.mt(c.mt(:,1) < c.mt_le_cen(1),:);
% 	phi = atan2(c.mt_le_cen(2) - mt(:,2),c.mt_le_cen(1) - mt(:,1)) * 360 / (2*pi);
% 	
% 	% Calculate the angle the stagnation point makes with the leading edge
% 	[phi,i] = unique(phi); mt = mt(i,:);
% 	dtdx = grad_mg(mt(:,1),mt(:,2));
% 	psi = atand(dtdx);
% 	
% 	% Correct the angle by the leading edge metal angle
% 	p.psi_stag = 90 - (c.chi_le - interp1(phi,psi,phi_stag,'pchip'));
% 	psi_stag = p.psi_stag;
% 	
% 	% iteration number
% 	ni = 1;
% 	beta=0.05;
% 	while abs(psi_stag-psi_target)>1.5
% 		[~,~]=system(['rm ' directory 'polarx.mises'],'-echo');
% 		if ni==1
% 			job.inlet_cam
% 			psi_stag
% 			psi_target
% 			dinlet_cam = (psi_stag-psi_target)*beta
% 			if abs(dinlet_cam) > 0.25
% 				dinlet_cam = 0.25*sign(dinlet_cam);
% 			end
% 			inlet_cam = job.inlet_cam + dinlet_cam;
% 			job.inlet_cam_0 = job.inlet_cam;
% 			job.psi_0 = psi_stag;
% 			job.inlet_cam=inlet_cam;
% 			
% 			[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle.py"','-echo');
% 			% change total camber to keep same exit angle
% 			job.tot_cam_pre = job.tot_cam;
% 			job.tot_cam = job.tot_cam - dinlet_cam;
% 			% change percentage ss camber to keep ss section the
% 			% same
% % % 			job.percentage_cam = double(double(job.percentage_cam) * double(job.tot_cam_pre)/double(job.tot_cam));
% 			[mt,c] = mises_create_mt(job,1);
% 			% Write section file
% 			save([directory 'section.mat'],'c','f');
% 			
% 		else
% 			%provide a limiter to increase in kappa_in
% 		    psi_stag
% 		    job.psi_0
% 			psi_target
%             dinlet_cam = (job.inlet_cam-job.inlet_cam_0)/(psi_stag-job.psi_0)*(psi_stag-psi_target)
% 			if abs(dinlet_cam) > 0.25
% 				dinlet_cam = 0.25*sign(dinlet_cam);
% 			end
% 			inlet_cam = job.inlet_cam - dinlet_cam;
% 			job.inlet_cam_0 = job.inlet_cam;
% 			job.psi_0 = psi_stag;
% 			job.inlet_cam = inlet_cam;
% 			job.tot_cam_pre = job.tot_cam;
% 			% 		% change total camber to keep same exit angle
% 			job.tot_cam = job.tot_cam - dinlet_cam;
% 			% change percentage ss camber to keep ss section the
% 			% same
% % 			job.percentage_cam = double(double(job.percentage_cam) * double(job.tot_cam_pre)/double(job.tot_cam));
% 			[mt,c] = mises_create_mt(job,1);
% 			% Write section file
% 			save([directory 'section.mat'],'c','f');
% 			[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle.py"','-echo');
% 		end
% 		
% 		ni = ni+1;
% 		
% 		
% 		
% 		% Read in flow file and re-run if not converged
% 		if exist([directory 'polarx.mises'],'file') ~= 0
% 			[Polarx, Ises] = mis_read_polarx('mises',directory);
% 			if isstruct(Polarx) ~= 1
% 				disp('*******Run again*********');
% 				while isstruct(Polarx) ~= 1 && n_type < 5
% 					n_type = n_type + 1;
% 					disp(['****Create new gridpar ' num2str(n_type) '****']);
% 					mis_gridpar(directory,n_type);
% 					[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle.py"','-echo');
% 					[Polarx, Ises] = mis_read_polarx('mises',directory);
% 				end
% 				n_type = 1; mis_gridpar(directory,n_type);
% 			end
% 		else
% 			disp('File Not Found')
% 			disp('*******Run again*********');
% 			Polarx = [];
% 			while isstruct(Polarx) ~= 1 && n_type < 5
% 				n_type = n_type + 1;
% 				disp(['****Create new gridpar ' num2str(n_type) '****']);
% 				mis_gridpar(directory,n_type);
% 				[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle.py"','-echo');
% 				[Polarx, Ises] = mis_read_polarx('mises',directory);
% 			end
% 			n_type = 1; mis_gridpar(directory,n_type);
% 			p = [];
% 			% 			return
% 		end
% 		
% 		% Read in the grid coodinates
% 		if exist([directory 'idat.mises_01'],'file') ~= 0
% 			Idat = mis_read_idat('mises_01',directory);
% 		else
% 			Idat = mis_read_idat('mises',directory);
% 		end
% 		
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
% 	end
% 	
% 	mises_directory = strrep(directory,'TURBOSTREAM','MISES');
% 	cd(mises_directory)
% 	% save new job
% 	save([mises_directory 'job.mat'],'job');
% end
% %% 
% 
% %% save new job.mat (will have changed due to incidence change)
% mises_directory = strrep(directory,'TURBOSTREAM','MISES');
% cd(mises_directory)
% % save new job
% save([mises_directory 'job.mat'],'job');
% 
% %% change ises back to normal smoothing
% % Write ises file
% fid = fopen([directory 'ises.mises'],'w');
% % fprintf(fid,'%s\n','1 2 5  !Global variables');
% % fprintf(fid,'%s\n','1 3 4  !Global constraints');
% % fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
% % fprintf(fid,'%s\n',' 3 4 6 9 17 !Global constraints');
% fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
% fprintf(fid,'%s\n',' 1 3 4 6 15 !Global constraints');
% 
% 
% fprintf(fid,'%6.5f %6.5f %6.5f %6.5f %6.5f ',[f.M f.Pin/f.Poin tand(f.Alpha) -0.66*m_chord f.Vtao1]);
% fprintf(fid,'%s\n','|Minl p1/p01 Sinl Xinl v1/ao1');
% % fprintf(fid,'0.00000 0.00000 0.00000 %6.5f |Mout p2/p01 Sout Xout\n',1.33*m_chord);
% fprintf(fid,'%6.5f %6.5f %6.5f %6.5f  ', [f.M2 f.Pout/f.Poin tand(f.Alpha2) 1.33*m_chord]);
% fprintf(fid,'%s\n','|Mout p2/p01 Sout Xout');
% fprintf(fid,'%s\n','0.00000     0.00000     1.39433     0.00000  | MFR  HWRATin GAMin DIFACTOR');
% fprintf(fid,'%i ',round(f.Re));
% fprintf(fid,'%s\n', '-4.00000 0.00000  0.26133 |Re -Turb');
% if is_turb == 1
%     fprintf(fid,'%s\n','0.001 0.01  |Xtr1 Xtr2');
% else
%     fprintf(fid,'%s\n','1.10 1.10  |Xtr1 Xtr2');
% end
% fprintf(fid,'%s\n','4     0.95000  1.00000  0.00000   | ISMOM  MCRIT  MUCON  PLOSSIN');
% fprintf(fid,'%s\n','0.00000  0.00000   | BVR1in  BVR2in');
% fclose(fid);
% 
% if run_single == 1
% 	[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingleHo-on.py"','-echo');
% else
% [~,~] = system('bash --login -c "python ~/bin/CalculateBladeMises_Ho-on.py"','-echo');	
% end
% 
% end


%% manufacture code (from Rig_260_Create)
% nb number of gum scanned blades
nb = size(d.thick,1);
n_type = 1;
% nb = 3;

% Run manufactured blade geometry
for n = 1:nb
	
	disp(['*****Running manufacturing number ' num2str(n) '****'])
	
	% load job.mat from intial MISES directory
	load([directory 'job.mat']);
	
	%% make directory
    directory_man = [directory(1:end-1) num2str(n) '/'];
%     directory_man = [directory(1:end-1) 'b' num2str(n) '/'];
	
	% Make directory if required
	if exist(directory_man,'dir') == 0
		mkdir(directory_man);
		[~,~]=system(['cp -r ' directory  '* ' directory_man]);
		[~,~]=system(['rm ' directory_man 'polarx.mises'],'-echo');
	else 
		continue;
	end
	
	% update directory
	job.man_directory = directory_man;
	
	
	load([directory_man 'section.mat'],'c','f');
	
	% Maximum MISES length
	m_chord = max(c.mt(:,1));
	
	%% change ises back to high smoothing
	% Write ises file
	fid = fopen([directory_man 'ises.mises'],'w');
	% fprintf(fid,'%s\n','1 2 5  !Global variables');
	% fprintf(fid,'%s\n','1 3 4  !Global constraints');
	% fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
	% fprintf(fid,'%s\n',' 3 4 6 9 17 !Global constraints');
	fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
	fprintf(fid,'%s\n',' 1 3 4 6 15 !Global constraints');
	
	
	fprintf(fid,'%6.5f %6.5f %6.5f %6.5f %6.5f ',[f.M f.Pin/f.Poin tand(f.Alpha) -0.66*m_chord f.Vtao1]);
	fprintf(fid,'%s\n','|Minl p1/p01 Sinl Xinl v1/ao1');
	% fprintf(fid,'0.00000 0.00000 0.00000 %6.5f |Mout p2/p01 Sout Xout\n',1.33*m_chord);
	fprintf(fid,'%6.5f %6.5f %6.5f %6.5f  ', [f.M2 f.Pout/f.Poin tand(f.Alpha2) 1.33*m_chord]);
	fprintf(fid,'%s\n','|Mout p2/p01 Sout Xout');
	fprintf(fid,'%s\n','0.00000     0.00000     1.39433     0.00000  | MFR  HWRATin GAMin DIFACTOR');
	fprintf(fid,'%i ',round(f.Re));
	fprintf(fid,'%s\n', '-4.00000 0.00000  0.26133 |Re -Turb');
	if is_turb == 1
		fprintf(fid,'%s\n','0.001 0.01  |Xtr1 Xtr2');
	else
		fprintf(fid,'%s\n','1.10 1.10  |Xtr1 Xtr2');
	end
	fprintf(fid,'%s\n','4     0.95000  -4.00000  0.00000   | ISMOM  MCRIT  MUCON  PLOSSIN');
	fprintf(fid,'%s\n','0.00000  0.00000   | BVR1in  BVR2in');
	fclose(fid);
	
	% re-define nominal thickness
	mf = t.mf;
	t = job.t;
	t.mf = mf;

	% extract mid_section thicknesses
	nj_mid = round(0.5*size(job.t.s_cl,2));
	if isfield(t,'thick_chord') == 1
		t_nom.thick_max= job.chord*t.thick_chord ;
	else
		t_nom.thick_max = t.thick_max(nj_mid);
	end
	t_nom.thick = interp1(t.s_cl(:,nj_mid),t.thick(:,nj_mid),t.mf{n,1}.s_cl,'pchip');

	
	% Record spanwise variations in thickness
	dev{1}.r_nondim = [0.1 0.5 0.9]; dev{1}.s = t.mf{n,1}.s_cl;
	for j = 1:size(d.thick,2)
		dev{1}.thick(:,j) = d.thick{n,j} * t_nom.thick_max + ...
			d.thick_max{n,j} * t_nom.thick;
		if isfield(d,'thick_sides_ss') == 1 && isfield(d,'thick_sides_ps') == 1
			dev{1}.thick_ss(:,j) = d.thick_sides_ss{n,j}* t_nom.thick_max;
			dev{1}.thick_ps(:,j) = d.thick_sides_ps{n,j}* t_nom.thick_max;
		end
		
		if isfield(d,'xy_le') == 1
			dev{1}.xy_le(:,j) = d.xy_le{n,j};
		end
	end
	
	% Smooth thickness distributions
	for j = 1:size(d.thick,2)
		
		dev{1}.thick(:,j) = smooth(dev{1}.s,dev{1}.thick(:,j));
		dev{1}.thick(1,j) = 0;
		dev{1}.thick(:,j) = smooth(dev{1}.thick(:,j),25);
		if isfield(d,'thick_sides_ss') == 1 && isfield(d,'thick_sides_ps') == 1
			dev{1}.thick_ss(:,j) = smooth(dev{1}.s,dev{1}.thick_ss(:,j));
			dev{1}.thick_ss(1,j) = 0;
			dev{1}.thick_ss(:,j) = smooth(dev{1}.thick_ss(:,j),25);
			dev{1}.thick_ps(:,j) = smooth(dev{1}.s,dev{1}.thick_ps(:,j));
			dev{1}.thick_ps(1,j) = 0;
			dev{1}.thick_ps(:,j) = smooth(dev{1}.thick_ps(:,j),25);
		end
	end

    %% pick which r_nondim to use
	section =  3; % is non-dimensional at 90% 
	%%
	man_dev.thick = dev{1}.thick(:,section);
	if isfield(d,'thick_sides_ss') == 1 && isfield(d,'thick_sides_ps') == 1
	man_dev.thick_ss = dev{1}.thick_ss(:,section);
	man_dev.thick_ps = dev{1}.thick_ps(:,section);
	end
	if isfield(d,'xy_le') == 1
		man_dev.xy_le = dev{1}.xy_le(:,section);
	end
	man_dev.s = dev{1}.s;
	job.dev = man_dev; % need to plot deviations later on
	
    % apply deviations due to chi_le;
% 	job.inlet_cam = job.inlet_cam + d.chi_le{n,section};
% 	job.tot_cam = job.tot_cam + d.chi_le{n,section}  - d.chi_te{n,section};
    
	% create new mt coordinates for mises
	write_file = 1;
	[mt,c] = mises_create_mt(job,write_file,man_dev);
	% Write section file
	save([directory_man 'section.mat'],'c','f');
	
    save([directory_man 'job.mat'],'job');
	
	cd(directory_man)
	
	[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle.py"','-echo');
	
	
	% Read in flow file and re-run if not converged
	if exist([directory_man 'polarx.mises'],'file') ~= 0
		[Polarx, Ises] = mis_read_polarx('mises',directory_man);
		if isstruct(Polarx) ~= 1
			disp('*******Run again*********');
			while isstruct(Polarx) ~= 1 && n_type < 5
				n_type = n_type + 1;
				disp(['****Create new gridpar ' num2str(n_type) '****']);
				mis_gridpar(directory_man,n_type);
				[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle.py"','-echo');
				[Polarx, Ises] = mis_read_polarx('mises',directory_man);
			end
			if n_type>5
				disp('*******Have to run manually*********');
				while isstruct(Polarx) ~= 1
					[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle_manual.py"','-echo');
					[Polarx, Ises] = mis_read_polarx('mises',directory_man);
				end
			end
			n_type = 1; mis_gridpar(directory_man,n_type);
		end
	else
		disp('File Not Found')
		disp('*******Run again*********');
		Polarx = [];
		while isstruct(Polarx) ~= 1 && n_type < 5
			n_type = n_type + 1;
			disp(['****Create new gridpar ' num2str(n_type) '****']);
			mis_gridpar(directory_man,n_type);
			[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle.py"','-echo');
			[Polarx, Ises] = mis_read_polarx('mises',directory_man);
		end
		if n_type==5
			disp('*******Have to run manually*********');
			while isstruct(Polarx) ~= 1
				[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle_manual.py"','-echo');
				[Polarx, Ises] = mis_read_polarx('mises',directory_man);
			end
		end
		n_type = 1; mis_gridpar(directory_man,n_type);
		p = [];
		% 			return
	end
	
	
	%% change ises back to normal smoothing
	% Write ises file
	fid = fopen([directory_man 'ises.mises'],'w');
	% fprintf(fid,'%s\n','1 2 5  !Global variables');
	% fprintf(fid,'%s\n','1 3 4  !Global constraints');
	% fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
	% fprintf(fid,'%s\n',' 3 4 6 9 17 !Global constraints');
	fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
	fprintf(fid,'%s\n',' 1 3 4 6 15 !Global constraints');
	
	
	fprintf(fid,'%6.5f %6.5f %6.5f %6.5f %6.5f ',[f.M f.Pin/f.Poin tand(f.Alpha) -0.66*m_chord f.Vtao1]);
	fprintf(fid,'%s\n','|Minl p1/p01 Sinl Xinl v1/ao1');
	% fprintf(fid,'0.00000 0.00000 0.00000 %6.5f |Mout p2/p01 Sout Xout\n',1.33*m_chord);
	fprintf(fid,'%6.5f %6.5f %6.5f %6.5f  ', [f.M2 f.Pout/f.Poin tand(f.Alpha2) 1.33*m_chord]);
	fprintf(fid,'%s\n','|Mout p2/p01 Sout Xout');
	fprintf(fid,'%s\n','0.00000     0.00000     1.39433     0.00000  | MFR  HWRATin GAMin DIFACTOR');
	fprintf(fid,'%i ',round(f.Re));
	fprintf(fid,'%s\n', '-4.00000 0.00000  0.26133 |Re -Turb');
	if is_turb == 1
		fprintf(fid,'%s\n','0.001 0.01  |Xtr1 Xtr2');
	else
		fprintf(fid,'%s\n','1.10 1.10  |Xtr1 Xtr2');
	end
	fprintf(fid,'%s\n','4     0.95000  1.00000  0.00000   | ISMOM  MCRIT  MUCON  PLOSSIN');
	fprintf(fid,'%s\n','0.00000  0.00000   | BVR1in  BVR2in');
	fclose(fid);
	
	%% run ises and polar
	if run_single == 1
	[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingleHo-on.py"','-echo');
	if exist('job','var') == 0
	[p,h] = mis_plot_section(directory_man,h,col,plot_stuff);
	else
		[p,h] = mis_plot_section(directory_man,h,col,plot_stuff,job);
	end
	else
		[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMises_Ho-on.py"','-echo');
		[Polar]=mis_read_polar('mises',directory_man);
		hold on;
		grid on
		plot(Polar.ainl-abs(f.Alpha),Polar.yp,'Color',col);
		ylabel(['Y_{p}'],'Rotation',0)
		xlabel(['\Delta \alpha'])
		p=[];
	end
	
end


end

