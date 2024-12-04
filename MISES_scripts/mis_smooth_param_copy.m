function [p,h] = mis_smooth_param(directory,h,col,plot_stuff,is_turb,job)
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

if exist('prior_solution','var') == 0
	prior_solution = 1;
end

if exist('job','var') == 1
% 	[c]=mises_output_xrrt(job);
% 	[f]=mises_output_finput(job);
% 	g=ts_read_hdf5([job.rjm_directory job.outname]);
% 	[~,~,delta] = ts_incidence_2D(g,0);
	psi_target = 0;
end

color = [0,0,1;1,0,0;0,1,0;0,0.5,0;1,0,1;0.5,0,0.5;0.5,0.5,0.5;0,0,0];

directory = strrep(directory,'TURBOSTREAM','MISES');

[p,h] = mis_plot_section(directory,[],[0,0,0],plot_stuff,job);
[p.Area1,p.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory,[],0);

%% create target
% put a limit on the value of A0 used to mormalise A10
if p.Area1<0.005
	A0 = 0.005;
else
	A0 = p.Area1;
end
if p.Area2<0.008
	A20 = 0.008;
else
	A20 = p.Area2;
end
A_TE01 = p.aHb_1;
A_TE02 = p.aHb_2;

if isfield(job,'Hk')==1
	H_TE = job.Hk;
else
	% 	H_TE = p.Hb_te;
% 	H_TE = 1.85
end
psi_target = 0;

% obtain alpha_target from loading coefficient
alpha_target = 41.5;
alpha_target = -atand(1/job.Vx_ratio*(-tand(job.alpha_rel_inlet)-job.load_coef/job.throttle_phi));

M_factor = [];

% pick s/c
% s_c = 1.2
% job.N_high = 2*pi*job.high_radius/(s_c*job.chord);

dHTE_ref=0.1;
dpsi_ref=5;
dalpha_ref=2;

% target = [-p.Area1/A0 ; -p.Area2/A20;  -p.aHb_1/A_TE01; -p.aHb_2/A_TE02; (H_TE-p.Hb_te)/dHTE_ref; (psi_target - p.psi_stag)/dpsi_ref]
% target = [-p.Area1/A0 ; -p.Area2/A20;  -p.aHb_1/A_TE01; -p.aHb_2/A_TE02; (H_TE-p.Hb_te)/dHTE_ref; (psi_target - p.psi_stag)/dpsi_ref; (alpha_target - p.Alpha)/dalpha_ref]
target = [-p.Area1/A0 ; -p.Area2/A20;  -p.aHb_1/A_TE01; -p.aHb_2/A_TE02; (psi_target - p.psi_stag)/dpsi_ref; (alpha_target - p.Alpha)/dalpha_ref]

if p.aHb_1<5*10^-4
	target(3) = 0;
end
if p.aHb_2<5*10^-4
	target(4) = 0;
end


%%

%% manufacture code (from Rig_260_Create)
% nb number of iterations
nb =1;
n_inner = 3;
% Run manufactured blade geometry
for n = 1:nb
	
	if n>8
		
		color(n,:) = [0,0,1];
	end
	
	
	disp(['*****Running global iteration number ' num2str(n) '****'])
	
	% load job.mat from intial MISES directory
	if n==1
% 		load([directory 'job.mat']);
	else
		load([directory_smooth 'job_final_' num2str(n-1) '.mat']);
	end
	
	%% make directory
% 	directory_smooth = [directory(1:end-1) num2str(n) '/'];
	directory_smooth = [directory(1:end-1) 'b' num2str(n) '/'];
	
	% Make directory if required
	if exist(directory_smooth,'dir') == 0
		mkdir(directory_smooth);
		[~,~]=system(['cp -r ' directory  '* ' directory_smooth]);
		[~,~]=system(['rm ' directory_smooth 'polarx.mises'],'-echo');
	else
		if n==1
		[p,h] = mis_plot_section(directory_smooth,[],[0,0,0],plot_stuff,job);
		else
			[p,h] = mis_plot_section(directory_smooth,h,color(n,:),plot_stuff,job);
		end
		[p.Area1,p.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
		[directory_smooth 'job_final_' num2str(n)  '.mat']
		load([directory_smooth 'job_final_' num2str(n)  '.mat']);
		target = job.p.target
		continue;
	end
	
	for n_i = 1:n_inner
		
		% update directory
		job.directory_smooth = directory_smooth;
		save([directory_smooth 'job.mat'],'job');
		
		load([directory_smooth 'section.mat'],'c','f');
% 		f.M = 0.7;
	    save(['section.mat'],'c','f');
		
		% Maximum MISES length
		m_chord = max(c.mt(:,1));
		
		%% change ises back to high smoothing
		% Write ises file
		fid = fopen([directory_smooth 'ises.mises'],'w');
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
		
		
		% create new mt coordinates for mises
		write_file = 1;
		[mt,c] = mises_create_mt(job,write_file);
		% Write section file
		save([directory_smooth 'section.mat'],'c','f');
		
		job.knot_coef = job.knot_coef + [0 0 0.1 0 0 0 0 0];
		save([directory_smooth 'job1.mat'],'job');
		[mt,c] = mises_create_mt(job,write_file);
		mis_run_inner(directory_smooth)
		[p1,h] = mis_plot_section(directory_smooth,h,[0,0,1],0,job);
		[p1.Area1,p1.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
		load([directory_smooth 'job.mat']);
		[~,~]=system(['rm ' directory_smooth 'polarx.mises']);
		[~,~]=system(['rm ' directory_smooth 'idat.mises*']);
		
		job.knot_coef = job.knot_coef + [0 0 0 0.1 0 0 0 0];
		save([directory_smooth 'job2.mat'],'job');
		[mt,c] = mises_create_mt(job,write_file);
		mis_run_inner(directory_smooth)
		[p2,h] = mis_plot_section(directory_smooth,h,[1,0,0],0,job);
		[p2.Area1,p2.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
		load([directory_smooth 'job.mat']);
		[~,~]=system(['rm ' directory_smooth 'polarx.mises']);
		[~,~]=system(['rm ' directory_smooth 'idat.mises*']);
		
		job.knot_coef = job.knot_coef + [0 0 0 0 0 0.15 0 0];
		save([directory_smooth 'job3.mat'],'job');
		[mt,c] = mises_create_mt(job,write_file);
		mis_run_inner(directory_smooth)
		[p3,h] = mis_plot_section(directory_smooth,h,[0,0.5,0],0,job);
		[p3.Area1,p3.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
		load([directory_smooth 'job.mat']);
		[~,~]=system(['rm ' directory_smooth 'polarx.mises']);
		[~,~]=system(['rm ' directory_smooth 'idat.mises*']);
		
% 		job.knot_coef = job.knot_coef + [0 0 0 0 0 0 0.15 0];
		job.knot_coef = job.knot_coef + [0 0 0 0 0.10 0 0 0];
		save([directory_smooth 'job4.mat'],'job');
		[mt,c] = mises_create_mt(job,write_file);
		mis_run_inner(directory_smooth)
		[p4,h] = mis_plot_section(directory_smooth,h,[0.5,0,0],0,job);
		[p4.Area1,p4.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
		load([directory_smooth 'job.mat']);
		[~,~]=system(['rm ' directory_smooth 'polarx.mises']);
		[~,~]=system(['rm ' directory_smooth 'idat.mises*']);
		
% 		job.N_high = job.N_high/1.05;
% % 		job.knot_coef = job.knot_coef + [0 0 0 0 0.10 0 0 0];
% 		save([directory_smooth 'job5.mat'],'job');
% 		[mt,c] = mises_create_mt(job,write_file);
% 		mis_run_inner(directory_smooth)
% 		[p5,h] = mis_plot_section(directory_smooth,h,[0,0,0.5],0,job);
% 		[p5.Area1,p5.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
% 		load([directory_smooth 'job.mat']);
% 		[~,~]=system(['rm ' directory_smooth 'polarx.mises']);
% 		[~,~]=system(['rm ' directory_smooth 'idat.mises*']);
		
		job.inlet_cam = job.inlet_cam+1.0;
		job.tot_cam = job.tot_cam + 1.0;
		save([directory_smooth 'job6.mat'],'job');
		[mt,c] = mises_create_mt(job,write_file);
		mis_run_inner(directory_smooth)
% 		[p6,h] = mis_plot_section(directory_smooth,h,[0.5,0.5,0.5],0,job);
% 		[p6.Area1,p6.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
		[p5,h] = mis_plot_section(directory_smooth,h,[0.5,0.5,0.5],0,job);
		[p5.Area1,p5.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
		load([directory_smooth 'job.mat']);
		[~,~]=system(['rm ' directory_smooth 'polarx.mises']);
		[~,~]=system(['rm ' directory_smooth 'idat.mises*']);
		
	    job.tot_cam = job.tot_cam+1.0;
		save([directory_smooth 'job7.mat'],'job');
		[mt,c] = mises_create_mt(job,write_file);
		mis_run_inner(directory_smooth)
		[p6,h] = mis_plot_section(directory_smooth,h,[0.5,0.5,0.5],0,job);
		[p6.Area1,p6.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
% 		[p7,h] = mis_plot_section(directory_smooth,h,[0.5,0.5,0.5],0,job);
% 		[p7.Area1,p7.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
		load([directory_smooth 'job.mat']);
		[~,~]=system(['rm ' directory_smooth 'polarx.mises']);
		[~,~]=system(['rm ' directory_smooth 'idat.mises*']);
		
% 		dR_dQ = [(p1.Area1-p.Area1)/A0 	(p2.Area1-p.Area1)/A0 (p3.Area1-p.Area1)/A0 (p4.Area1-p.Area1)/A0 (p5.Area1-p.Area1)/A0 (p6.Area1-p.Area1)/A0; ...
% 			(p1.Area2-p.Area2)/A20 	(p2.Area2-p.Area2)/A20 (p3.Area2-p.Area2)/A20 (p4.Area2-p.Area2)/A20 (p5.Area2-p.Area2)/A20 (p6.Area2-p.Area2)/A20; ...
% 			(p1.aHb_1-p.aHb_1)/A_TE01 	(p2.aHb_1-p.aHb_1)/A_TE01 (p3.aHb_1-p.aHb_1)/A_TE01 (p4.aHb_1-p.aHb_1)/A_TE01 (p5.aHb_1-p.aHb_1)/A_TE01 (p6.aHb_1-p.aHb_1)/A_TE01; ...
% 			(p1.aHb_2-p.aHb_2)/A_TE02 	(p2.aHb_2-p.aHb_2)/A_TE02 (p3.aHb_2-p.aHb_2)/A_TE02 (p4.aHb_2-p.aHb_2)/A_TE02 (p5.aHb_2-p.aHb_2)/A_TE02 (p6.aHb_2-p.aHb_2)/A_TE02; ...
% 			(p1.Hb_te-p.Hb_te)/dHTE_ref 	(p2.Hb_te-p.Hb_te)/dHTE_ref (p3.Hb_te-p.Hb_te)/dHTE_ref (p4.Hb_te-p.Hb_te)/dHTE_ref (p5.Hb_te-p.Hb_te)/dHTE_ref (p6.Hb_te-p.Hb_te)/dHTE_ref; ...
% 			(p1.psi_stag-p.psi_stag)/dpsi_ref 	(p2.psi_stag-p.psi_stag)/dpsi_ref (p3.psi_stag-p.psi_stag)/dpsi_ref (p4.psi_stag-p.psi_stag)/dpsi_ref (p5.psi_stag-p.psi_stag)/dpsi_ref (p6.psi_stag-p.psi_stag)/dpsi_ref;
% 			]
		
% 		dR_dQ = [(p1.Area1-p.Area1)/A0 	(p2.Area1-p.Area1)/A0 (p3.Area1-p.Area1)/A0 (p4.Area1-p.Area1)/A0 (p5.Area1-p.Area1)/A0 (p6.Area1-p.Area1)/A0 (p7.Area1-p.Area1)/A0; ...
% 			(p1.Area2-p.Area2)/A20 	(p2.Area2-p.Area2)/A20 (p3.Area2-p.Area2)/A20 (p4.Area2-p.Area2)/A20 (p5.Area2-p.Area2)/A20 (p6.Area2-p.Area2)/A20 (p7.Area2-p.Area2)/A20; ...
% 			(p1.aHb_1-p.aHb_1)/A_TE01 	(p2.aHb_1-p.aHb_1)/A_TE01 (p3.aHb_1-p.aHb_1)/A_TE01 (p4.aHb_1-p.aHb_1)/A_TE01 (p5.aHb_1-p.aHb_1)/A_TE01 (p6.aHb_1-p.aHb_1)/A_TE01 (p7.aHb_1-p.aHb_1)/A_TE01; ...
% 			(p1.aHb_2-p.aHb_2)/A_TE02 	(p2.aHb_2-p.aHb_2)/A_TE02 (p3.aHb_2-p.aHb_2)/A_TE02 (p4.aHb_2-p.aHb_2)/A_TE02 (p5.aHb_2-p.aHb_2)/A_TE02 (p6.aHb_2-p.aHb_2)/A_TE02 (p7.aHb_2-p.aHb_2)/A_TE02;...
% 			(p1.Hb_te-p.Hb_te)/dHTE_ref 	(p2.Hb_te-p.Hb_te)/dHTE_ref (p3.Hb_te-p.Hb_te)/dHTE_ref (p4.Hb_te-p.Hb_te)/dHTE_ref (p5.Hb_te-p.Hb_te)/dHTE_ref (p6.Hb_te-p.Hb_te)/dHTE_ref (p7.Hb_te-p.Hb_te)/dHTE_ref; ...
% 			(p1.psi_stag-p.psi_stag)/dpsi_ref 	(p2.psi_stag-p.psi_stag)/dpsi_ref (p3.psi_stag-p.psi_stag)/dpsi_ref (p4.psi_stag-p.psi_stag)/dpsi_ref (p5.psi_stag-p.psi_stag)/dpsi_ref (p6.psi_stag-p.psi_stag)/dpsi_ref (p7.psi_stag-p.psi_stag)/dpsi_ref;
% 			(p1.Alpha-p.Alpha)/dalpha_ref 	(p2.Alpha-p.Alpha)/dalpha_ref (p3.Alpha-p.Alpha)/dalpha_ref (p4.Alpha-p.Alpha)/dalpha_ref (p5.Alpha-p.Alpha)/dalpha_ref (p6.Alpha-p.Alpha)/dalpha_ref (p7.Alpha-p.Alpha)/dalpha_ref;
% 			]

		dR_dQ = [(p1.Area1-p.Area1)/A0 	(p2.Area1-p.Area1)/A0 (p3.Area1-p.Area1)/A0 (p4.Area1-p.Area1)/A0 (p5.Area1-p.Area1)/A0 (p6.Area1-p.Area1)/A0; ...
			(p1.Area2-p.Area2)/A20 	(p2.Area2-p.Area2)/A20 (p3.Area2-p.Area2)/A20 (p4.Area2-p.Area2)/A20 (p5.Area2-p.Area2)/A20 (p6.Area2-p.Area2)/A20; ...
			(p1.aHb_1-p.aHb_1)/A_TE01 	(p2.aHb_1-p.aHb_1)/A_TE01 (p3.aHb_1-p.aHb_1)/A_TE01 (p4.aHb_1-p.aHb_1)/A_TE01 (p5.aHb_1-p.aHb_1)/A_TE01 (p6.aHb_1-p.aHb_1)/A_TE01; ...
			(p1.aHb_2-p.aHb_2)/A_TE02 	(p2.aHb_2-p.aHb_2)/A_TE02 (p3.aHb_2-p.aHb_2)/A_TE02 (p4.aHb_2-p.aHb_2)/A_TE02 (p5.aHb_2-p.aHb_2)/A_TE02 (p6.aHb_2-p.aHb_2)/A_TE02;...
			(p1.psi_stag-p.psi_stag)/dpsi_ref 	(p2.psi_stag-p.psi_stag)/dpsi_ref (p3.psi_stag-p.psi_stag)/dpsi_ref (p4.psi_stag-p.psi_stag)/dpsi_ref (p5.psi_stag-p.psi_stag)/dpsi_ref (p6.psi_stag-p.psi_stag)/dpsi_ref;
			(p1.Alpha-p.Alpha)/dalpha_ref 	(p2.Alpha-p.Alpha)/dalpha_ref (p3.Alpha-p.Alpha)/dalpha_ref (p4.Alpha-p.Alpha)/dalpha_ref (p5.Alpha-p.Alpha)/dalpha_ref (p6.Alpha-p.Alpha)/dalpha_ref;
			]
		
		cond=rcond(dR_dQ)
		dQ_dR = inv(dR_dQ);
		
		% calculate change in dQ
		dQ = dQ_dR * target
		
% 		dQ = [0.15; 0.15; 0.15; 0.15; 0.05*2*pi*job.high_radius/job.N_high/job.chord; 1.0].*dQ;
% 		dQ = [0.10; 0.10; 0.15; 0.10; 0.05*2*pi*job.high_radius/job.N_high/job.chord; 1.0; 1.0].*dQ;
	    dQ = [0.10; 0.10; 0.15; 0.10; 1.0; 1.0].*dQ;
		job.dQ = dQ;
		
		% provide limits to change in dQ
		if abs(job.dQ(1,1)) > 0.10
			job.dQ(1,1) = 0.10*sign(job.dQ(1,1));
		end
		if abs(job.dQ(2,1)) > 0.10
			job.dQ(2,1) = 0.10*sign(job.dQ(2,1));
		end
		if abs(job.dQ(3,1)) > 0.15
			job.dQ(3,1) = 0.15*sign(job.dQ(3,1));
		end
		if abs(job.dQ(4,1)) > 0.10
			job.dQ(4,1) = 0.10*sign(job.dQ(4,1));
		end
		% 		if abs(job.dQ(5,1)) > 0.05*2*pi*job.high_radius/job.N_high/job.chord
		% 			job.dQ(5,1) = 0.05*2*pi*job.high_radius/job.N_high/job.chord*sign(job.dQ(5,1));
		% 		end
		if abs(job.dQ(5,1)) > 1
			job.dQ(5,1) = 1*sign(job.dQ(5,1));
		end
		if abs(job.dQ(6,1)) > 1
			job.dQ(6,1) = 1*sign(job.dQ(6,1));
		end
		
		dQ = job.dQ
		
		dQ(isnan(job.dQ)==1)=0;
		job.dQ = dQ;
		
		% run final blade
		% 		job.knot_coef = job.knot_coef + [ 0 0 job.dQ(1,1) job.dQ(2,1) 0 job.dQ(3,1) job.dQ(4,1) 0];
		job.knot_coef = job.knot_coef + [ 0 0 job.dQ(1,1) job.dQ(2,1) job.dQ(4,1) job.dQ(3,1) 0 0];
% 		N_factor = job.dQ(5,1)/(2*pi*job.high_radius/job.N_high/job.chord)
% 		job.N_high = job.N_high/(1+N_factor);
% 		job.inlet_cam = job.inlet_cam + job.dQ(6,1);
% 		jb.tot_cam = job.tot_cam + job.dQ(7,1);
	    job.inlet_cam = job.inlet_cam + job.dQ(5,1);
	    job.tot_cam = job.tot_cam + job.dQ(6,1);
		
		% for controlled diffusion blades
		if job.knot_coef(3)>1.1
			job.knot_coef(3) = 1;
		end
		
		save([directory_smooth 'job_final_' num2str(n) '.mat'],'job');
		[mt,c] = mises_create_mt(job,write_file);
		mis_run_inner(directory_smooth)
		[p,h] = mis_plot_section(directory_smooth,h,[1,0,0],plot_stuff,job);
		[p.Area1,p.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
		
		%% change ises back to normal smoothing
		% Write ises file
		fid = fopen([directory_smooth 'ises.mises'],'w');
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
		
		mis_run_inner(directory_smooth)
% 		[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingleHo-on.py"','-echo');
		[p,h] = mis_plot_section(directory_smooth,h,color(n_i,:),plot_stuff,job);
		[p.Area1,p.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
		
		%% create new target
		% put a limit on the value of A0 used to mormalise A10
		if p.Area1<0.005
			A0 = 0.005;
		else
			A0 = p.Area1;
		end
		if p.Area2<0.008
			A20 = 0.008;
		else
			A20 = p.Area2;
		end
		A_TE01 = p.aHb_1;
		A_TE02 = p.aHb_2;
		
% 		target = [-p.Area1/A0 ; -p.Area2/A20;  -p.aHb_1/A_TE01; -p.aHb_2/A_TE02; (H_TE-p.Hb_te)/dHTE_ref; (psi_target - p.psi_stag)/dpsi_ref]
% 		target = [-p.Area1/A0 ; -p.Area2/A20;  -p.aHb_1/A_TE01; -p.aHb_2/A_TE02; (H_TE-p.Hb_te)/dHTE_ref; (psi_target - p.psi_stag)/dpsi_ref; (alpha_target - p.Alpha)/dalpha_ref]
		target = [-p.Area1/A0 ; -p.Area2/A20;  -p.aHb_1/A_TE01; -p.aHb_2/A_TE02; (psi_target - p.psi_stag)/dpsi_ref; (alpha_target - p.Alpha)/dalpha_ref]
		
		if p.aHb_1<3*10^-3
			target(3) = 0;
		end
		if p.aHb_2<3*10^-3
			target(4) = 0;
		end
		
		p.target = target;
		job.p  = p;
		save([directory_smooth 'job_final_' num2str(n) '.mat'],'job');
		
		%%
		
		% check whether solution has converged
% 		if abs(p.Area1)<0.005 && abs(p.Area2)<0.0075 && abs(p.aHb_1)<2*10^-3 && abs(p.aHb_2)<2*10^-3 && abs(p.Hb_te-H_TE)<0.02 && abs(p.psi_stag)<2
% 		if abs(p.Area1)<0.005 && abs(p.Area2)<0.0075 && abs(p.aHb_1)<2*10^-3 && abs(p.aHb_2)<2*10^-3 && abs(p.Hb_te-H_TE)<0.02 && abs(p.psi_stag)<2	&& abs(p.Alpha-alpha_target)<1.0
		if abs(p.Area1)<0.005 && abs(p.Area2)<0.0075 && abs(p.aHb_1)<2*10^-3 && abs(p.aHb_2)<2*10^-3 && abs(p.psi_stag)<2 && abs(p.Alpha-alpha_target)<1.0	
			disp('*****Reached smoothness convergence*****')
			break
		end
		
		%% correct for local incidence if off by a lot
		if abs(p.psi_stag) > 3.5
			disp('****Have to restagger****');
			psi_stag = p.psi_stag
			% iteration number
			ni = 1;
		    n_psi_max = 6;
			beta=0.05;
			n_type = 1;
			n_type_lim = 6; 
			while abs(psi_stag-psi_target)>2 && ni<n_psi_max
				[~,~]=system(['rm ' directory_smooth 'polarx.mises'],'-echo');
				[~,~]=system(['rm ' directory_smooth 'idat.mises*']);
				if ni==1
					job.inlet_cam
					psi_stag
					psi_target
					dinlet_cam = (psi_stag-psi_target)*beta
					if abs(dinlet_cam) > 1.0
						dinlet_cam = 1.0*sign(dinlet_cam);
					end
					inlet_cam = job.inlet_cam + dinlet_cam;
					job.inlet_cam_0 = job.inlet_cam;
					% 			job.psi_0 = delta.psi;
					job.psi_0 = psi_stag;
					job.inlet_cam=inlet_cam;
					
					% change total camber to keep same exit angle
					job.tot_cam_pre = job.tot_cam;
					job.tot_cam = job.tot_cam - dinlet_cam;
					% change percentage ss camber to keep ss section the
					% same
					% 			job.percentage_cam = double(double(job.percentage_cam) * double(job.tot_cam_pre)/double(job.tot_cam));
					[mt,c] = mises_create_mt(job,1);
					% Write section file
					save([directory_smooth 'section.mat'],'c','f');
					[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle.py"','-echo');
				else
					%provide a limiter to increase in kappa_in
					psi_target
					psi_stag
					job.psi_0
					% 			dinlet_cam = (job.inlet_cam-job.inlet_cam_0)/(delta.psi-job.psi_0)*(delta.psi-job.psi_target)
					dinlet_cam = (job.inlet_cam-job.inlet_cam_0)/(psi_stag-job.psi_0)*(psi_stag-psi_target)
					if abs(dinlet_cam) > 1.0
						dinlet_cam = 1.0*sign(dinlet_cam);
					end
					inlet_cam = job.inlet_cam - dinlet_cam;
					job.inlet_cam_0 = job.inlet_cam;
					% 			job.psi_0 = delta.psi;
					job.psi_0 = psi_stag;
					job.inlet_cam = inlet_cam;
					job.tot_cam_pre = job.tot_cam;
					% 		% change total camber to keep same exit angle
					job.tot_cam = job.tot_cam - dinlet_cam;
					% change percentage ss camber to keep ss section the
					% same
					% 			job.percentage_cam = double(double(job.percentage_cam) * double(job.tot_cam_pre)/double(job.tot_cam));
					[mt,c] = mises_create_mt(job,1);
					% Write section file
					save([directory_smooth 'section.mat'],'c','f');
					[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle.py"','-echo');
					
				end
				
				ni = ni+1;
				
				% Read in flow file and re-run if not converged
				if exist([directory_smooth 'polarx.mises'],'file') ~= 0
					[Polarx, Ises] = mis_read_polarx('mises',directory_smooth);
					if isstruct(Polarx) ~= 1
						disp('*******Run again*********');
						while isstruct(Polarx) ~= 1 && n_type < n_type_lim
							n_type = n_type + 1;
							disp(['****Create new gridpar ' num2str(n_type) '****']);
							mis_gridpar(directory_smooth,n_type);
							[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle.py"','-echo');
							[Polarx, Ises] = mis_read_polarx('mises',directory_smooth);
						end
						if n_type>4
							disp('*******Have to run manually*********');
							while isstruct(Polarx) ~= 1
								[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle_manual.py"','-echo');
								[Polarx, Ises] = mis_read_polarx('mises',directory);
							end
						end
						n_type = 1; mis_gridpar(directory,n_type);
					end
				else
					disp('File Not Found')
					disp('*******Run again*********');
					Polarx = [];
					while isstruct(Polarx) ~= 1 && n_type < n_type_lim
						n_type = n_type + 1;
						disp(['****Create new gridpar ' num2str(n_type) '****']);
						mis_gridpar(directory_smooth,n_type);
						[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle.py"','-echo');
						[Polarx, Ises] = mis_read_polarx('mises',directory_smooth);
					end
					if n_type>4
						disp('*******Have to run manually*********');
						while isstruct(Polarx) ~= 1
							[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMisesSingle_manual.py"','-echo');
							[Polarx, Ises] = mis_read_polarx('mises',directory);
						end
					end
					n_type = 1; mis_gridpar(directory,n_type);
					p = [];
					% 			return
				end
				
				[p,h] = mis_plot_section(directory_smooth,h,[1,0,0],0,job);
				psi_stag = p.psi_stag;
			end
			
			%%
			[p,h] = mis_plot_section(directory_smooth,h,[1,0,0],plot_stuff,job);
			[p.Area1,p.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
			
			%% create new target
			% put a limit on the value of A0 used to mormalise A10
			if p.Area1<0.005
				A0 = 0.005;
			else
				A0 = p.Area1;
			end
			if p.Area2<0.008
				A20 = 0.008;
			else
				A20 = p.Area2;
			end
			A_TE01 = p.aHb_1;
			A_TE02 = p.aHb_2;
			
% 			target = [-p.Area1/A0 ; -p.Area2/A20;  -p.aHb_1/A_TE01; -p.aHb_2/A_TE02; (H_TE-p.Hb_te)/dHTE_ref; (psi_target - p.psi_stag)/dpsi_ref]
% 			target = [-p.Area1/A0 ; -p.Area2/A20;  -p.aHb_1/A_TE01; -p.aHb_2/A_TE02; (H_TE-p.Hb_te)/dHTE_ref; (psi_target - p.psi_stag)/dpsi_ref; (alpha_target - p.Alpha)/dalpha_ref]
			target = [-p.Area1/A0 ; -p.Area2/A20;  -p.aHb_1/A_TE01; -p.aHb_2/A_TE02; (psi_target - p.psi_stag)/dpsi_ref; (alpha_target - p.Alpha)/dalpha_ref]
			
			if p.aHb_1<3*10^-3
				target(3) = 0;
			end
			if p.aHb_2<3*10^-3
				target(4) = 0;
			end
			
			p.target = target;
			job.p  = p;
			save([directory_smooth 'job_final_' num2str(n) '.mat'],'job');
			
		end
		
		if abs(p.Area1)<0.005 && abs(p.Area2)<0.0075 && abs(p.aHb_1)<2*10^-3 && abs(p.aHb_2)<2*10^-3 && abs(p.psi_stag)<2 && abs(p.Alpha-alpha_target)<1.0
			disp('*****Reached smoothness convergence*****')
			break
		end
		
	end
	
	
end

