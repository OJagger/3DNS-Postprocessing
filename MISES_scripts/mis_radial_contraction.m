function [rad_contr,AtA1,AtA1_roV,rad_contr_o,throat_min,xrt_throat] = mis_radial_contraction(directory,plot_stuff,bl_inc)
% rad_contr: average radial contraction
% rad_contr_o: average radial contraction along the throat

if exist('plot_stuff','var')==0 || isempty(plot_stuff) == 1
	plot_stuff = 0;
end


directory = strrep(directory,'TURBOSTREAM','MISES');
directory = strrep(directory,'/home/dl467/','/mnt/Disk1/dl467/');
mis_directory = directory;

cd(mis_directory)
load([mis_directory 'job.mat']);
load([mis_directory 'section.mat']);


% extract stagnation line from mises
% Read in flow file

if exist([mis_directory 'idat.mises_01'],'file') ~= 0
    Idat = mis_read_idat('mises_01',mis_directory);
else
    Idat = mis_read_idat('mises',mis_directory);
end

% read stream file
[stream] = mis_read_stream('mises',mis_directory);

if isfield(job,'r_le') == 1
    pitch = 2*pi*0.5*(c.r_le+c.r_te)/c.N;
else
    pitch = 2*pi*(job.r_mid)/c.N;
end
if isfield(job,'tchord')==0
    job.tchord = job.chord;
end
s_c = pitch/job.tchord;

%% if you want to calculate with boundary layer %%
if exist('bl_inc','var')==0 || isempty(bl_inc) == 1 || bl_inc == 0
	bl=[];
elseif bl_inc == 1
bl = [Idat.x(:,1).*c.r_le Idat.y(:,1).*c.r_le+pitch];
x = bl(:,1);
y = bl(:,2);
% construct section
ote = 1; dev =[];
cam_spline = 1;
% construct section
%% this is legacy
if isfield(c,'xs_sb')==0
    [xrt,sl_cx,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
    xrt_ss = mca.xrt_ss; xrt_ps = mca.xrt_ps;
    %%
else
    [xrt,xrt_cam,xrt_ps,xrt_ss]=bl_construct_section(job,0,ote,dev,cam_spline);
end

y = y(x>min(xrt_ss(:,1)/pitch));
x = x(x>min(xrt_ss(:,1)/pitch));
bl = [x,y];
end
%%
	

% xrt_stag_line = [Idat.x(1:Idat.ninl(2),end)*c.r_le , Idat.y(1:Idat.ninl(2),end)*c.r_le];
xrt_stag_line = [Idat.x(:,end) Idat.y(:,end)].*c.r_le;
plot_throat = 0;
% pitch_ref = job.chord;
pitch_ref = pitch;
%% use normals on midpassage streamline to get throat
mode = directory;
%%

% plot it from geometry
alpha_in = abs(job.alpha_rel_inlet);

if plot_stuff == 1
h1 = figure(); grid on; hold on;
h2 = figure(); grid on; hold on;
figure(); grid on; hold on;
plot_throat = 1;
end
% construct section
ote = 1; dev =[];
cam_spline = 1;
% construct section
%% this is legacy
if isfield(c,'xs_sb')==0
    [xrt,sl_cx,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
    xrt_ss = mca.xrt_ss; xrt_ps = mca.xrt_ps;
    %%
else
    [xrt,xrt_cam,xrt_ps,xrt_ss]=bl_construct_section(job,0,ote,dev,cam_spline);
end

if isfield(job,'alpha_rel_inlet')==1
    [~,throat_min,xrt_throat,xrt_planes] = ts_throat_calc(xrt_ss,xrt_ps,pitch,bl,plot_throat,pitch_ref,xrt_stag_line,alpha_in,mode);
else
    [~,throat_min,xrt_throat,~] = ts_throat_calc(xrt_ss,xrt_ps,pitch,bl,plot_throat,pitch_ref,xrt_stag_line);
end

for m=1:size(xrt_planes,2)
	if isempty(xrt_planes{m})==0
	x_cx(m) = xrt_planes{m}.x(1);
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

[x_inter,rt_inter]=intersections(s_stream_x*c.r_le/pitch,s_stream_rt*c.r_le/pitch+1,xrt_throat(:,1),xrt_throat(:,2));
% calculate intersection point on suction surface
[x_inter2,rt_inter2]=intersections(s_stream_ss_x*c.r_le/pitch,s_stream_ss_rt*c.r_le/pitch+1,xrt_throat(:,1),xrt_throat(:,2));



[x_inter,~] = unique(x_inter);
[~,i_passage]=min(abs(s_stream_x*c.r_le/pitch-x_inter));

if isempty(x_inter2) == 0
    [x_inter2,~] = unique(x_inter2);
    if size(x_inter2,1) ~= 1
        x_inter2 = x_inter2(1);
    end
    [~,i_passage2]=min(abs(s_stream_ss_x*c.r_le/pitch-x_inter2));
end


[~,i_le] = min(abs(s_stream_ss_x));
stream_contr_ss = interp1(stream.x,stream.b,s_stream_ss_x,'linear');
[~,i_te] = min(stream_contr_ss);


s_stream_throat(1) = s_stream_x(i_passage)*c.r_le;
s_stream_throat(2) = s_stream_rt(i_passage)*c.r_le;
s_stream_throat(3) = (s_stream_ss(i_passage)-s_stream_ss(i_le))/(s_stream_ss(i_te)-s_stream_ss(i_le));
s_stream_throat(3) = 0.5*(s_stream_throat(3)+1);

% calculate midpassage streamline contraction
rad_contr = stream_contr(i_passage);

% calculate radial contraction along throat
rad_contr_o = interp1(stream.x*c.r_le/pitch,stream.b,xrt_throat(:,1));

Idat.xt = xrt_throat(:,1);
Idat.yt = xrt_throat(:,2);
Idat.et = -1./[ xrt_throat(end,1)-xrt_throat(1,1) , xrt_throat(end,2)-xrt_throat(1,2)];
Idat.et = Idat.et./sqrt(Idat.et(:,1).^2+Idat.et(:,2).^2);
Idat.etn(:,1) = -Idat.et(:,2);
Idat.etn(:,2) = Idat.et(:,1);
Idat.xt_cell = 0.5*(xrt_throat(1:end-1,1)+xrt_throat(2:end,1));
Idat.yt_cell = 0.5*(xrt_throat(1:end-1,2)+xrt_throat(2:end,2));

% calculate normal area to streamline of each cell
% take radial contraction into account
Idat.Ax_ratio = interp1(stream.x*c.r_le/pitch,stream.b,Idat.xt_cell);
dAx = diff(Idat.yt).*Idat.Ax_ratio;
dAt = diff(Idat.xt).*Idat.Ax_ratio;
dA_t = sqrt(dAx.^2+dAt.^2);
At = sum(dA_t);
A1 = cosd(Idat.binl);
AtA1 = At./A1;

% calculate area ratio based on roV
% get important primary variables along throat
varnames = {'Vx','Vt','P','ro','q','Po','dPo_Po','M'};
Idat.ro = Idat.ro(1:end-1,1:end-1);
for j = 1:length(varnames)
    throat.([varnames{j}]) = griddata(Idat.x_cell.*c.r_le./pitch_ref,Idat.y_cell.*c.r_le./pitch_ref+1,Idat.(varnames{j}),Idat.xt_cell,Idat.yt_cell);
end
throat.V = sqrt(throat.Vx.^2+throat.Vt.^2);
es = [ throat.Vx./throat.V , throat.Vt./throat.V];
etn = repmat(Idat.etn,size(es,1),1);
dm = throat.V.*throat.ro.*dAx;
throat.roV_mass = sum(throat.V.*sum(etn.*es,2).*throat.ro.*dm)./sum(dm);
throat.roV_area = sum(throat.V.*sum(etn.*es,2).*throat.ro.*dA_t)./sum(dA_t);
roV1_mass=sum(Idat.q1.*Idat.ro(1,:).*Idat.dm_tot)./(sum(Idat.dm_tot));
AtA1_roV = roV1_mass./throat.roV_mass;

if plot_stuff == 1
plot(Idat.x(1:end-1,:)*c.r_le/pitch_ref,Idat.y(1:end-1,:)*c.r_le/pitch_ref+1,'-k');
plot(s_stream_x*c.r_le/pitch_ref,s_stream_rt*c.r_le/pitch_ref+1,'-r');

plot(s_stream_ss_x(i_le)*c.r_le/pitch_ref,s_stream_ss_rt(i_le)*c.r_le/pitch_ref+1,'xr');
plot(s_stream_ss_x(i_te)*c.r_le/pitch_ref,s_stream_ss_rt(i_te)*c.r_le/pitch_ref+1,'xr');
plot(s_stream_x(i_passage)*c.r_le/pitch_ref,s_stream_rt(i_passage)*c.r_le/pitch_ref+1,'xr');
if isempty(x_inter2) == 0
plot(s_stream_ss_x(i_passage2)*c.r_le/pitch_ref,s_stream_ss_rt(i_passage2)*c.r_le/pitch_ref+1,'xg');
end

plot(xrt_stag_line(:,1)/pitch_ref,xrt_stag_line(:,2)/pitch_ref,'-b')
plot(xrt_throat(:,1),xrt_throat(:,2),'-k')
plot(s_stream_x(i_passage)*c.r_le/pitch,s_stream_rt(i_passage)*c.r_le/pitch+1,'xg')

figure(h1); grid on; hold on;
plot(s_stream_x,stream_contr)
plot(s_stream_x(i_passage),stream_contr(i_passage),'xr')

figure(h2); grid on; hold on;
plot(roV1_mass./(throat.V.*sum(etn.*es,2).*throat.ro),cumsum(dm./sum(dm)));
ylabel('Non-dimensional mass-flow'); xlabel('At/A1 local')
end

end