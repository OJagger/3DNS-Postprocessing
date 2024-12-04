function [p,h,job] = mis_smooth_param2(directory,oscosa_tar,col,plot_stuff,mode,job,alpha_target)
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

% mode = 1 - change s/c
% mode = 2 - change percentage_camber

% Prepare defaults
if exist('col','var') == 0 || isempty(col) == 1
	col = [0 0 0];
end
if exist('plot_stuff','var') == 0
	plot_stuff = 1;
end
% if exist('is_turb','var') == 0
% 	is_turb = 1;
% end

is_turb =1;

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

controlled_diffusion = 0;

%% create target
% put a limit on the value of A0 used to mormalise A10
if abs(p.Area1)<0.005
	A0 = 0.005;
else
	A0 = p.Area1;
end
if abs(p.Area2)<0.008
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
if exist('alpha_target','var') == 0 || isempty(alpha_target) == 1
alpha_target = 41.5;
alpha_target = abs(-atand(1/job.Vx_ratio*(-tand(job.alpha_rel_inlet)-job.load_coef/job.throttle_phi)))
end

M_factor = [];

[xrt,sl_cx,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
% calculate value of o_scosa
[~,o_s]=ts_throat_calc(mca.xrt_ss,mca.xrt_ps,2*pi*job.high_radius/job.N_high,[],0);
p.oscosa=o_s/cosd(job.alpha_rel_inlet)
if isfield(job,'Mt_des')==1
	[o_scosa_des]=Mt_oscosa(job.Mt_des,job)
else
	o_scosa_des = p.oscosa;
end

if isempty(oscosa_tar) == 1
if abs(p.oscosa-o_scosa_des)<0.05
	oscosa_target = p.oscosa
% 	oscosa_target = 1.084;
else
	oscosa_target = o_scosa_des
% 	oscosa_target = 1.084;
end
else
	oscosa_target = oscosa_tar
end

p.oscosa_target = oscosa_target;

% pick s/c
% s_c = 1.2
% job.N_high = 2*pi*job.high_radius/(s_c*job.chord);

dHTE_ref=0.1;
dpsi_ref=5;
dalpha_ref=2;
doscosa_ref=0.025*oscosa_target;

% target = [-p.Area1/A0 ; -p.Area2/A20;  -p.aHb_1/A_TE01; -p.aHb_2/A_TE02; (H_TE-p.Hb_te)/dHTE_ref; (psi_target - p.psi_stag)/dpsi_ref]
% target = [-p.Area1/A0 ; -p.Area2/A20;  -p.aHb_1/A_TE01; -p.aHb_2/A_TE02; (H_TE-p.Hb_te)/dHTE_ref; (psi_target - p.psi_stag)/dpsi_ref; (alpha_target - p.Alpha)/dalpha_ref]
target = [-p.Area1/A0 ; -p.Area2/A20;  -p.aHb_1/A_TE01; -p.aHb_2/A_TE02; (psi_target - p.psi_stag)/dpsi_ref; (alpha_target - p.Alpha)/dalpha_ref; (oscosa_target-p.oscosa)/doscosa_ref];

if abs(p.aHb_1)<1*10^-3
	target(3) = 0;
end
if abs(p.aHb_2)<1*10^-3
	target(4) = 0;
end

if abs(p.Area1)<0.005
	target(1) = 0;
end
if abs(p.Area2)<0.01
	target(2) = 0;
end

target

frac_chord_org = job.pitch_chord;
percentage_cam_org = job.percentage_cam;
knot_coef_org = job.knot_coef(5);

%%
flag = 0;

%% manufacture code (from Rig_260_Create)
% nb number of iterations
nb =4;
n_inner = 3;
% Run manufactured blade geometry
for n = 1:nb
	
	if(flag==1)
		flag
		break
	end
	
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
	directory_smooth = [directory(1:end-1) 'c' num2str(n) '/'];
	
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
		[xrt,sl_cx,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
		% calculate value of o_scosa
		[~,o_s]=ts_throat_calc(mca.xrt_ss,mca.xrt_ps,2*pi*job.high_radius/job.N_high,[],0);
		p.oscosa=o_s/cosd(job.alpha_rel_inlet)
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
		fprintf(fid,'%s\n','4     0.65000  -4.00000  0.00000   | ISMOM  MCRIT  MUCON  PLOSSIN');
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
		[xrt,sl_cx,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
		% calculate value of o_scosa
		[~,o_s]=ts_throat_calc(mca.xrt_ss,mca.xrt_ps,2*pi*job.high_radius/job.N_high,[],0);
		p1.oscosa=o_s/cosd(job.alpha_rel_inlet);
		load([directory_smooth 'job.mat']);
		[~,~]=system(['rm ' directory_smooth 'polarx.mises']);
		[~,~]=system(['rm ' directory_smooth 'idat.mises*']);
		
		job.knot_coef = job.knot_coef + [0 0 0 0.1 0 0 0 0];
		save([directory_smooth 'job2.mat'],'job');
		[mt,c] = mises_create_mt(job,write_file);
		mis_run_inner(directory_smooth)
		[p2,h] = mis_plot_section(directory_smooth,h,[1,0,0],0,job);
		[p2.Area1,p2.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
		[xrt,sl_cx,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
		% calculate value of o_scosa
		[~,o_s]=ts_throat_calc(mca.xrt_ss,mca.xrt_ps,2*pi*job.high_radius/job.N_high,[],0);
		p2.oscosa=o_s/cosd(job.alpha_rel_inlet);
		load([directory_smooth 'job.mat']);
		[~,~]=system(['rm ' directory_smooth 'polarx.mises']);
		[~,~]=system(['rm ' directory_smooth 'idat.mises*']);
		
		job.knot_coef = job.knot_coef + [0 0 0 0 0 0.10 0 0];
		save([directory_smooth 'job3.mat'],'job');
		[mt,c] = mises_create_mt(job,write_file);
		mis_run_inner(directory_smooth)
		[p3,h] = mis_plot_section(directory_smooth,h,[0,0.5,0],0,job);
		[p3.Area1,p3.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
		[xrt,sl_cx,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
		% calculate value of o_scosa
		[~,o_s]=ts_throat_calc(mca.xrt_ss,mca.xrt_ps,2*pi*job.high_radius/job.N_high,[],0);
		p3.oscosa=o_s/cosd(job.alpha_rel_inlet);
		load([directory_smooth 'job.mat']);
		[~,~]=system(['rm ' directory_smooth 'polarx.mises']);
		[~,~]=system(['rm ' directory_smooth 'idat.mises*']);
		
		job.knot_coef = job.knot_coef + [0 0 0 0 0 0 0.10 0];
		% 		job.knot_coef = job.knot_coef + [0 0 0 0 0.10 0 0 0];
		save([directory_smooth 'job4.mat'],'job');
		[mt,c] = mises_create_mt(job,write_file);
		mis_run_inner(directory_smooth)
		[p4,h] = mis_plot_section(directory_smooth,h,[0.5,0,0],0,job);
		[p4.Area1,p4.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
		[xrt,sl_cx,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
		% calculate value of o_scosa
		[~,o_s]=ts_throat_calc(mca.xrt_ss,mca.xrt_ps,2*pi*job.high_radius/job.N_high,[],0);
		p4.oscosa=o_s/cosd(job.alpha_rel_inlet);
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
		
		job.inlet_cam = job.inlet_cam+0.5;
		job.tot_cam = job.tot_cam + 0.5;
		save([directory_smooth 'job6.mat'],'job');
		[mt,c] = mises_create_mt(job,write_file);
		mis_run_inner(directory_smooth)
		% 		[p6,h] = mis_plot_section(directory_smooth,h,[0.5,0.5,0.5],0,job);
		% 		[p6.Area1,p6.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
		[p5,h] = mis_plot_section(directory_smooth,h,[0.5,0.5,0.5],0,job);
		[p5.Area1,p5.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
		[xrt,sl_cx,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
		% calculate value of o_scosa
		[~,o_s]=ts_throat_calc(mca.xrt_ss,mca.xrt_ps,2*pi*job.high_radius/job.N_high,[],0);
		p5.oscosa=o_s/cosd(job.alpha_rel_inlet);
		load([directory_smooth 'job.mat']);
		[~,~]=system(['rm ' directory_smooth 'polarx.mises']);
		[~,~]=system(['rm ' directory_smooth 'idat.mises*']);
		
		job.tot_cam_temp = job.tot_cam+1.0;
		job.percentage_cam_temp = job.percentage_cam * job.tot_cam/job.tot_cam_temp;
		job.tot_cam = job.tot_cam_temp;
		save([directory_smooth 'job7.mat'],'job');
		[mt,c] = mises_create_mt(job,write_file);
		mis_run_inner(directory_smooth)
		[p6,h] = mis_plot_section(directory_smooth,h,[0.5,0.5,0.5],0,job);
		[p6.Area1,p6.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
		[xrt,sl_cx,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
		% calculate value of o_scosa
		[~,o_s]=ts_throat_calc(mca.xrt_ss,mca.xrt_ps,2*pi*job.high_radius/job.N_high,[],0);
		p6.oscosa=o_s/cosd(job.alpha_rel_inlet);
		% 		[p7,h] = mis_plot_section(directory_smooth,h,[0.5,0.5,0.5],0,job);
		% 		[p7.Area1,p7.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
		load([directory_smooth 'job.mat']);
		[~,~]=system(['rm ' directory_smooth 'polarx.mises']);
		[~,~]=system(['rm ' directory_smooth 'idat.mises*']);
		
		% 		job.knot_coef = job.knot_coef + [0 0 0 0 0 0 0.10 0];
		
		if mode == 1
			job.pitch_chord = job.pitch_chord + 0.05;
		elseif mode == 2
			job.knot_coef = job.knot_coef + [0 0 0 0 0.05 0 0 0];
%             job.percentage_cam = job.percentage_cam + 0.05;
		end
		%
		save([directory_smooth 'job8.mat'],'job');
		[mt,c] = mises_create_mt(job,write_file);
		mis_run_inner(directory_smooth)
		[p7,h] = mis_plot_section(directory_smooth,h,[0.5,0,0],0,job);
		[p7.Area1,p7.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
		[xrt,sl_cx,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
		% calculate value of o_scosa
		[~,o_s]=ts_throat_calc(mca.xrt_ss,mca.xrt_ps,2*pi*job.high_radius/job.N_high,[],0);
		p7.oscosa=o_s/cosd(job.alpha_rel_inlet);
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
		
		dR_dQ = [(p1.Area1-p.Area1)/A0 	(p2.Area1-p.Area1)/A0 (p3.Area1-p.Area1)/A0 (p4.Area1-p.Area1)/A0 (p5.Area1-p.Area1)/A0 (p6.Area1-p.Area1)/A0 (p7.Area1-p.Area1)/A0;...
			(p1.Area2-p.Area2)/A20 	(p2.Area2-p.Area2)/A20 (p3.Area2-p.Area2)/A20 (p4.Area2-p.Area2)/A20 (p5.Area2-p.Area2)/A20 (p6.Area2-p.Area2)/A20 (p7.Area2-p.Area2)/A20; ...
			(p1.aHb_1-p.aHb_1)/A_TE01 	(p2.aHb_1-p.aHb_1)/A_TE01 (p3.aHb_1-p.aHb_1)/A_TE01 (p4.aHb_1-p.aHb_1)/A_TE01 (p5.aHb_1-p.aHb_1)/A_TE01 (p6.aHb_1-p.aHb_1)/A_TE01 (p7.aHb_1-p.aHb_1)/A_TE01; ...
			(p1.aHb_2-p.aHb_2)/A_TE02 	(p2.aHb_2-p.aHb_2)/A_TE02 (p3.aHb_2-p.aHb_2)/A_TE02 (p4.aHb_2-p.aHb_2)/A_TE02 (p5.aHb_2-p.aHb_2)/A_TE02 (p6.aHb_2-p.aHb_2)/A_TE02 (p7.aHb_2-p.aHb_2)/A_TE02;...
			(p1.psi_stag-p.psi_stag)/dpsi_ref 	(p2.psi_stag-p.psi_stag)/dpsi_ref (p3.psi_stag-p.psi_stag)/dpsi_ref (p4.psi_stag-p.psi_stag)/dpsi_ref (p5.psi_stag-p.psi_stag)/dpsi_ref (p6.psi_stag-p.psi_stag)/dpsi_ref (p7.psi_stag-p.psi_stag)/dpsi_ref;
			(p1.Alpha-p.Alpha)/dalpha_ref 	(p2.Alpha-p.Alpha)/dalpha_ref (p3.Alpha-p.Alpha)/dalpha_ref (p4.Alpha-p.Alpha)/dalpha_ref (p5.Alpha-p.Alpha)/dalpha_ref (p6.Alpha-p.Alpha)/dalpha_ref (p7.Alpha-p.Alpha)/dalpha_ref; ...
			(p1.oscosa-p.oscosa)/doscosa_ref 	(p2.oscosa-p.oscosa)/doscosa_ref (p3.oscosa-p.oscosa)/doscosa_ref (p4.oscosa-p.oscosa)/doscosa_ref (p5.oscosa-p.oscosa)/doscosa_ref (p6.oscosa-p.oscosa)/doscosa_ref (p7.oscosa-p.oscosa)/doscosa_ref; ...
			]
		
		cond=rcond(dR_dQ)
		dQ_dR = inv(dR_dQ);
		
		% calculate change in dQ
		dQ = dQ_dR * target
		
		% 		dQ = [0.15; 0.15; 0.15; 0.15; 0.05*2*pi*job.high_radius/job.N_high/job.chord; 1.0].*dQ;
		% 		dQ = [0.10; 0.10; 0.15; 0.10; 0.05*2*pi*job.high_radius/job.N_high/job.chord; 1.0; 1.0].*dQ;
		dQ = [0.10; 0.10; 0.10; 0.10; 0.5; 1.0; 0.05].*dQ;
		job.dQ = dQ;
		
		% provide limits to change in dQ
		if abs(job.dQ(1,1)) > 0.10
			job.dQ(1,1) = 0.10*sign(job.dQ(1,1));
		end
		if abs(job.dQ(2,1)) > 0.10
			job.dQ(2,1) = 0.10*sign(job.dQ(2,1));
		end
		if abs(job.dQ(3,1)) > 0.10
			job.dQ(3,1) = 0.10*sign(job.dQ(3,1));
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
		if abs(job.dQ(7,1)) > 0.05
            job.dQ(7,1) = 0.05*sign(job.dQ(7,1));
		end
		
		dQ = job.dQ
		
		dQ(isnan(job.dQ)==1)=0;
		job.dQ = dQ;
		
		if mode == 2
		% limit change in percentage cam
		if abs(job.knot_coef(5) - knot_coef_org)>0.1
			job.dQ(7,1) = 0;
		end
		end
		if abs(job.percentage_cam - percentage_cam_org)>0.1
			job.dQ(7,1) = 0;
		end

		
		% run final blade
		if mode == 1
		job.knot_coef = job.knot_coef + [ 0 0 job.dQ(1,1) job.dQ(2,1) 0 job.dQ(3,1) job.dQ(4,1) 0];
		elseif mode == 2 
        job.knot_coef = job.knot_coef + [ 0 0 job.dQ(1,1) job.dQ(2,1) job.dQ(7,1) job.dQ(3,1) job.dQ(4,1) 0];
%           job.knot_coef = job.knot_coef + [ 0 0 job.dQ(1,1) job.dQ(2,1) 0 job.dQ(3,1) job.dQ(4,1) 0];
% 		  job.percentage_cam = job.percentage_cam + job.dQ(7,1);
		end
		job.knot_coef = abs(job.knot_coef);
		% 		job.knot_coef = job.knot_coef + [ 0 0 job.dQ(1,1) job.dQ(2,1) job.dQ(4,1) job.dQ(3,1) 0 0];
		% 		N_factor = job.dQ(5,1)/(2*pi*job.high_radius/job.N_high/job.chord)
		% 		job.N_high = job.N_high/(1+N_factor);
		% 		job.inlet_cam = job.inlet_cam + job.dQ(6,1);
		% 		jb.tot_cam = job.tot_cam + job.dQ(7,1);
		job.inlet_cam = job.inlet_cam + job.dQ(5,1);
		job.tot_cam_temp = job.tot_cam + job.dQ(6,1);
		job.percentage_cam = job.percentage_cam * job.tot_cam/job.tot_cam_temp;
		job.tot_cam = job.tot_cam_temp;
		if mode == 1
		job.pitch_chord_temp = job.pitch_chord + job.dQ(7,1);
		if abs(job.pitch_chord_temp-frac_chord_org)<0.10;
		job.pitch_chord = job.pitch_chord_temp;
		end
		end
		
		
		% check if camber line is too inflectional
		if job.knot_coef(4)/job.knot_coef(3)>1.3
			beta = 1.15;
			knot_coef_temp = job.knot_coef(3);
			job.knot_coef(3) = job.knot_coef(4)/beta;
			job.knot_coef(4) = knot_coef_temp*beta;
			disp('****Camber is becoming too inflectional***')
		elseif job.knot_coef(4)/job.knot_coef(3)<1/1.3
			beta = 1.15;
			knot_coef_temp = job.knot_coef(3);
			job.knot_coef(3) = job.knot_coef(4)*beta;
			job.knot_coef(4) = knot_coef_temp/beta;
			disp('****Camber is becoming too inflectional***')
		end
		if job.knot_coef(5)/job.knot_coef(4)>1.40 || job.knot_coef(5)/job.knot_coef(4)<1/1.40
			beta = 1.2;
			knot_coef_temp = job.knot_coef(5);
			job.knot_coef(5) = beta*job.knot_coef(4);
			job.knot_coef(4) = knot_coef_temp/beta;
			disp('****Camber is becoming too inflectional 2***')
		elseif job.knot_coef(5)/job.knot_coef(4)<1/1.40
			beta = 1.2;
			knot_coef_temp = job.knot_coef(5);
			job.knot_coef(5) = job.knot_coef(4)/beta;
			job.knot_coef(4) = knot_coef_temp*beta;
			disp('****Camber is becoming too inflectional 2***')
		end
		if job.knot_coef(7)/job.knot_coef(6)>1.45
			beta = 1.25;
			knot_coef_temp = job.knot_coef(7);
			job.knot_coef(7) = beta*job.knot_coef(6);
			job.knot_coef(6) = knot_coef_temp/beta;
			disp('****Camber is becoming too inflectional 3***')
		elseif job.knot_coef(7)/job.knot_coef(6)<1/1.40
			beta = 1.25;
			knot_coef_temp = job.knot_coef(7);
			job.knot_coef(7) = job.knot_coef(6)/beta;
			job.knot_coef(6) = knot_coef_temp*beta;
			disp('****Camber is becoming too inflectional 3***')
		end
		
		
		% check if inverse camber exists and remove it
		[xrt,~,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
		if min(mca.cam) < 0
			while min(mca.cam) < 0
				if mca.lamda1 ~= 0
					job.knot_coef = job.knot_coef + [ 0 0 -0.02 -0.01 0 0 0 0];
				else
					job.knot_coef = job.knot_coef + [ 0 0 0 -0.02 -0.01 0 0 0];
				end
				disp('****Inverse camber exists (remove it)****')
				[xrt,~,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
			end
		end
		
		% for controlled diffusion blades
		if controlled_diffusion == 1
			if job.knot_coef(3)>1.1
				job.knot_coef(3) = 1;
			end
		end
		
		save([directory_smooth 'job_final_' num2str(n) '.mat'],'job');
		[mt,c] = mises_create_mt(job,write_file);
		[p,h,job] = mis_run_section(directory_smooth,[],[0,0,1],0,1,1,job);
% 		mis_run_inner(directory_smooth)
		[p,h] = mis_plot_section(directory_smooth,h,[1,0,0],plot_stuff,job);
		[p.Area1,p.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
	    [xrt,sl_cx,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
		% calculate value of o_scosa
		[~,o_s]=ts_throat_calc(mca.xrt_ss,mca.xrt_ps,2*pi*job.high_radius/job.N_high,[],0);
		p.oscosa=o_s/cosd(job.alpha_rel_inlet);
		
		%% create new target
		% put a limit on the value of A0 used to mormalise A10
		if abs(p.Area1)<0.005
			A0 = 0.005;
		else
			A0 = p.Area1;
		end
		if abs(p.Area2)<0.008
			A20 = 0.008;
		else
			A20 = p.Area2;
		end
		A_TE01 = p.aHb_1;
		A_TE02 = p.aHb_2;
		
		% 		target = [-p.Area1/A0 ; -p.Area2/A20;  -p.aHb_1/A_TE01; -p.aHb_2/A_TE02; (H_TE-p.Hb_te)/dHTE_ref; (psi_target - p.psi_stag)/dpsi_ref]
		% 		target = [-p.Area1/A0 ; -p.Area2/A20;  -p.aHb_1/A_TE01; -p.aHb_2/A_TE02; (H_TE-p.Hb_te)/dHTE_ref; (psi_target - p.psi_stag)/dpsi_ref; (alpha_target - p.Alpha)/dalpha_ref]
% 		target = [-p.Area1/A0 ; -p.Area2/A20;  -p.aHb_1/A_TE01; -p.aHb_2/A_TE02; (psi_target - p.psi_stag)/dpsi_ref; (alpha_target - p.Alpha)/dalpha_ref]
		target = [-p.Area1/A0 ; -p.Area2/A20;  -p.aHb_1/A_TE01; -p.aHb_2/A_TE02; (psi_target - p.psi_stag)/dpsi_ref; (alpha_target - p.Alpha)/dalpha_ref; (oscosa_target-p.oscosa)/doscosa_ref];
		
		
		if abs(p.aHb_1)<5*10^-3
			target(3) = 0;
		end
		if abs(p.aHb_2)<5*10^-3
			target(4) = 0;
		end
		
		if abs(p.Area1)<0.005
			target(1) = 0;
		end
		if abs(p.Area2)<0.01
			target(2) = 0;
		end
		
		target
		p.target = target;
		job.p  = p;
		save([directory_smooth 'job_final_' num2str(n) '.mat'],'job');
		

		[p,h] = mis_plot_section(directory_smooth,h,[1,0,0],plot_stuff,job);
		[p.Area1,p.Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory_smooth,M_factor,0);
		[xrt,sl_cx,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
		% calculate value of o_scosa
		[~,o_s]=ts_throat_calc(mca.xrt_ss,mca.xrt_ps,2*pi*job.high_radius/job.N_high,[],0);
		p.oscosa=o_s/cosd(job.alpha_rel_inlet);
		%% create new target
		% put a limit on the value of A0 used to mormalise A10
		if abs(p.Area1)<0.005
			A0 = 0.005;
		else
			A0 = p.Area1;
		end
		if abs(p.Area2)<0.008
			A20 = 0.008;
		else
			A20 = p.Area2;
		end
		A_TE01 = p.aHb_1;
		A_TE02 = p.aHb_2;
		
		% 			target = [-p.Area1/A0 ; -p.Area2/A20;  -p.aHb_1/A_TE01; -p.aHb_2/A_TE02; (H_TE-p.Hb_te)/dHTE_ref; (psi_target - p.psi_stag)/dpsi_ref]
		% 			target = [-p.Area1/A0 ; -p.Area2/A20;  -p.aHb_1/A_TE01; -p.aHb_2/A_TE02; (H_TE-p.Hb_te)/dHTE_ref; (psi_target - p.psi_stag)/dpsi_ref; (alpha_target - p.Alpha)/dalpha_ref]
% 		target = [-p.Area1/A0 ; -p.Area2/A20;  -p.aHb_1/A_TE01; -p.aHb_2/A_TE02; (psi_target - p.psi_stag)/dpsi_ref; (alpha_target - p.Alpha)/dalpha_ref]
		target = [-p.Area1/A0 ; -p.Area2/A20;  -p.aHb_1/A_TE01; -p.aHb_2/A_TE02; (psi_target - p.psi_stag)/dpsi_ref; (alpha_target - p.Alpha)/dalpha_ref; (oscosa_target-p.oscosa)/doscosa_ref]
		
		
		if abs(p.aHb_1)<5*10^-3
			target(3) = 0;
		end
		if abs(p.aHb_2)<5*10^-3
			target(4) = 0;
		end
		
		if abs(p.Area1)<0.005
			target(1) = 0;
		end
		if abs(p.Area2)<0.01
			target(2) = 0;
		end
		
		target
		p.target = target;
		p.alpha_target = alpha_target;
		p.oscosa_target = oscosa_target;
		job.p  = p;
		save([directory_smooth 'job_final_' num2str(n) '.mat'],'job');
		
% 		if abs(p.Area1)<0.005 && abs(p.Area2)<0.02 && abs(p.aHb_1)<5*10^-3 && abs(p.aHb_2)<5*10^-3 && abs(p.psi_stag)<2 && abs(p.Alpha-alpha_target)<1.0 && abs(p.oscosa-oscosa_target)<0.0025
		% stricter restriction in alpha_target
        if abs(p.Area1)<0.005 && abs(p.Area2)<0.02 && abs(p.aHb_1)<5*10^-3 && abs(p.aHb_2)<5*10^-3 && abs(p.psi_stag)<2 && abs(p.Alpha-alpha_target)<0.15 && abs(p.oscosa-oscosa_target)<0.0025	
            disp('*****Reached smoothness convergence*****')
			flag = 1;
			break
		end
		
	end
	
	if(flag==1)
		flag
		break
	end
	
	
end

