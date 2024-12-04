% determine tthroat Mach number from mises simulation
function [Idat] = mis_throat_mach(directory,bl_inc,plot_stuff)

if exist('plot_stuff','var')==0 || isempty(plot_stuff) == 1
	plot_stuff = 0;
end


directory = strrep(directory,'TURBOSTREAM','MISES');
mis_directory = directory;

cd(mis_directory)
load([mis_directory 'job.mat']);
load([mis_directory 'section.mat']);


% Read in flow file
if exist([directory 'polarx.mises'],'file') ~= 0
    [Polarx, Ises] = mis_read_polarx('mises',directory);
else
    disp('File Not Found')
    p = [];
    return
end


% extract stagnation line from mises
% Read in flow file
if exist([mis_directory 'idat.mises_01'],'file') ~= 0
    Idat = mis_read_idat('mises_01',mis_directory);
else
    Idat = mis_read_idat('mises',mis_directory);
end

% read stream file
[stream] = mis_read_stream('mises',mis_directory);

%% legacy
if isfield(job,'r_le')==0
   job.r_le = job.r_mid; job.r_te = job.r_mid; 
end
%%

%% legacy
% get number of blades
if isfield(job,'N_high')==1
   pitch = 2*pi*(job.high_radius)/job.N_high;
   job.N = 1./(pitch./(pi*(job.r_le+job.r_te)));
end
%%
pitch = 2*pi*0.5*(job.r_le+job.r_te)/job.N;

if isfield(job,'tchord')==0
   job.tchord = job.chord; 
end
s_c = pitch/job.tchord;

%% if you want to calculate with boundary layer %%
if exist('bl_inc','var')==0 || isempty(bl_inc) == 1 || (bl_inc) == 0 
	bl=[];
else
bl = [Idat.x(:,1).*c.r_le Idat.y(:,1).*c.r_le];
x = bl(:,1);
y = bl(:,2);
% construct section WITH trailing edge to reparameterise
ote = 0; dev =[]; 
cam_spline = 1;
if isfield(job,'chi_le')==1
    [xrt,xrt_cam,xrt_ps,xrt_ss]=bl_construct_section(job,0,ote,dev,cam_spline);
    mca.xrt_ss = xrt_ss;
    mca.xrt_ps = xrt_ps;
else
    [xrt,sl_cx,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
end
y = y(x>min(mca.xrt_ss(:,1)))+pitch;
x = x(x>min(mca.xrt_ss(:,1)));
bl = [x,y];
end
%%
	

% xrt_stag_line = [Idat.x(1:Idat.ninl(2),end)*c.r_le , Idat.y(1:Idat.ninl(2),end)*c.r_le];
xrt_stag_line = [Idat.x(:,end) Idat.y(:,end)].*c.r_le;
plot_throat = 0;
% pitch_ref = job.chord;
pitch_ref = pitch;
mode = 1;

% plot it from geometry
alpha_in = abs(job.alpha_rel_inlet);

if plot_stuff == 1
figure; grid on; hold on;
contourf(Idat.x_cell.*c.r_le./pitch_ref,Idat.y_cell.*c.r_le./pitch_ref+1,Idat.M); colorbar;
plot_throat = 1;
end

% construct section WITH trailing edge to reparameterise
ote = 0; dev =[]; plot_stuff = 0;
cam_spline = 1;
if isfield(job,'chi_le')==1
    [xrt,xrt_cam,xrt_ps,xrt_ss]=bl_construct_section(job,plot_stuff,ote,dev,cam_spline);
    mca.xrt_ss = xrt_ss;
    mca.xrt_ps = xrt_ps;
else
    [xrt,sl_cx,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
end

if isfield(job,'alpha_rel_inlet')==1
    [throat,~,xrt_throat,xrt_throat_cov]=ts_throat_calc(mca.xrt_ss,mca.xrt_ps,pitch,bl,plot_throat,pitch_ref,xrt_stag_line,alpha_in,mode);
else
    [throat,~,xrt_throat,xrt_throat_cov]=ts_throat_calc(mca.xrt_ss,mca.xrt_ps,pitch,bl,plot_throat,pitch_ref,xrt_stag_line);
end

for m=1:size(xrt_throat_cov,2)
	if isempty(xrt_throat_cov{m})==0
	x_cx(m) = xrt_throat_cov{m}.x(1);
	else
		x_cx(m) = NaN;
	end
end


% plot it from mises for each streamtube
% Calculate blade normals to suction surface
for j = 1:size(Idat.x,2)-1
lower_streamline = [Idat.x(:,j) Idat.y(:,j)];
upper_streamline = [Idat.x(:,j+1) Idat.y(:,j+1)]; 
% streamwise distance
s_stream(:,j) = [0 ; cumsum(sum(diff(lower_streamline,1,1).^2,2).^0.5)];
s_stream(:,j) = (s_stream(:,j));
end

% calculate mid stream
s_stream_mid = s_stream(1:end-1,round(size(s_stream,2)/2));
s_stream_ss = s_stream(1:end-1,1);
s_stream_x = Idat.x(1:end-1,round(size(s_stream,2)/2));
s_stream_rt = Idat.y(1:end-1,round(size(s_stream,2)/2));
s_stream_ss_x = Idat.x(1:end-1,1);
s_stream_ss_rt= Idat.y(1:end-1,1);
stream_contr = interp1(stream.x,stream.b,s_stream_x,'linear');
 
[x_inter,rt_inter]=intersections(s_stream_x*c.r_le/pitch_ref,s_stream_rt*c.r_le/pitch_ref+1,xrt_throat(:,1),xrt_throat(:,2));
% calculate intersection point on suction surface
[x_inter2,rt_inter2]=intersections(s_stream_ss_x*c.r_le/pitch_ref,s_stream_ss_rt*c.r_le/pitch_ref+1,xrt_throat(:,1),xrt_throat(:,2));


[x_inter,~] = unique(x_inter);
[~,i_passage]=min(abs(s_stream_x*c.r_le/pitch_ref-x_inter));

if isempty(x_inter2) == 0
    [x_inter2,~] = unique(x_inter2);
    [~,i_passage2]=min(abs(s_stream_ss_x*c.r_le/pitch_ref-x_inter2(1)));
end

% calculate mass flow along throat

Idat.Vfrac_xy(:,:,1) = Idat.Vfrac.*Idat.es(:,:,1);
Idat.Vfrac_xy(:,:,2) = Idat.Vfrac.*Idat.es(:,:,2);
Idat.roVx=Idat.ro(1:end-1,1:end-1)./Idat.ro1.*Idat.Vfrac_xy(:,:,1)./(sum(Idat.es(1,:,1))./size(Idat.es,2));
for l = 1:size(xrt_throat,1)-1
Idat.xrt_cell(l,:) = 0.5*(xrt_throat(l,:)+xrt_throat(l+1,:));
end
Idat.s_t = [0 ; cumsum(sqrt(sum(diff(Idat.xrt_cell).^2,2)))];
s_t = Idat.s_t./Idat.s_t(end);
Idat.roVx_throat=griddata(Idat.x_cell.*c.r_le./pitch_ref,Idat.y_cell.*c.r_le./pitch_ref + 1,Idat.roVx,Idat.xrt_cell(:,1),Idat.xrt_cell(:,2));
Idat.Mt = griddata(Idat.x_cell.*c.r_le./pitch_ref,Idat.y_cell.*c.r_le./pitch_ref + 1,Idat.M,Idat.xrt_cell(:,1),Idat.xrt_cell(:,2));
Idat.dm_t = Idat.roVx_throat.*diff(xrt_throat(:,2));
Idat.Mtavg = sum(Idat.Mt.*Idat.dm_t)./(sum(Idat.dm_t));
if plot_stuff == 1
figure; plot(Idat.Mt,s_t)
end

end