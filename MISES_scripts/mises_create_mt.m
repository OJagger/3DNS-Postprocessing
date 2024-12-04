% create mt and parameterise blade c from job.mat file 
% option to write out blade file

function [mt,c] = mises_create_mt(job,write_file)


if exist('write_file','var')==0
	write_file = 1; % write blade.mises file by default
end

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

%% legacy
if isfield(job,'tchord')==0
    job.tchord = job.chord;
end
%%

% copy over thickness parameters from job
thick_names = {'rad_le','s_thick_max','rad_thick_max','wedge_te','thick_te','thick_max','tchord'};
for v = 1:length(thick_names)
    c.(thick_names{v}) = job.(thick_names{v});
end

pitch = 2*pi*0.5*(job.r_le+job.r_te)/job.N;
disp(['Pitch-to-chord is ' num2str(pitch/job.tchord) ])

% construct section WITH trailing edge to reparameterise
ote = 0; dev =[]; plot_stuff = 0;
cam_spline = 1;
[xrt,xrt_cam,xrt_ps,xrt_ss]=bl_construct_section(job,plot_stuff,ote,dev,cam_spline);


% reparameterise blade
% plot_stuff = 0; fast_mode = 0;
% [c,h] = bl_parameterise_section(xrt,fast_mode,plot_stuff,pitch);

% construct section WITHOUT trailing edge to input to mises
ote = 1; dev =[]; plot_stuff = 0;
cam_spline = 1; 
ni = 151; % number of points need to be reduced to run in mises
[xrt,xrt_cam,xrt_ps,xrt_ss]=bl_construct_section(job,plot_stuff,ote,dev,cam_spline,ni);
x = xrt(:,1); rt = xrt(:,2);
c.xrt_cam = xrt_cam;
c.s_cl = [0 ; cumsum(sum(diff(xrt_cam,1,1).^2,2).^0.5,1)]; 

bl = []; pitch_ref = pitch; plot_throat = 0;
[~,o_s,xrt_throat]=ts_throat_calc(xrt_ss,xrt_ps,pitch,bl,plot_throat,pitch_ref);
x_throat=xrt_throat(1,1).*pitch_ref;
s_cam = [0;cumsum((diff(xrt_cam(:,1),1).^2 + diff(xrt_cam(:,2),1).^2).^0.5)];
s_throat = interp1(xrt_cam(:,1),s_cam,x_throat);
c.ss_sb = s_throat./s_cam(end);
c.xs_sb = x_throat./(xrt_cam(end,1)-xrt_cam(1,1));

% Check blade nose always points down and inlet flow angle is positive
[~,i_min] = min(x);
if rt(i_min) > rt(1); rt = -rt; end; 


% Calculate MISES coordinates
r = 0.5*(job.r_le+job.r_te);
c.r_le = job.r_le; c.r_te = job.r_te;
c.N = job.N;
mt = [x ./ r rt ./ r]; 

% Maximum MISES length
m_chord = max(mt(:,1));
tchord = sum(sum(diff(xrt_cam,1,1).^2,2).^0.5,1);

% Record leading edge point for plotting later
c.mt_le = xrt_cam(1,:) / c.r_le;

% Interpolate the leading edge centre for plotting later
xrt_le_cen = [interp1(c.s_cl,c.xrt_cam(:,1),0.03,'pchip') interp1(c.s_cl,c.xrt_cam(:,2),0.03,'pchip')];
c.mt_le_cen = xrt_le_cen / interp1(c.xrt_cam([1 end],1),[c.r_le c.r_te],xrt_le_cen(1));

% save chi_le and chi_te coordinates
c.chi_le = atand(diff(xrt_cam(1:2,2),1,1) ./ diff(xrt_cam(1:2,1),1,1));
c.chi_te = atand(diff(xrt_cam(end-1:end,2),1,1) ./ diff(xrt_cam(end-1:end,1),1,1));

%% Save remaining coordinates to section data structure
c.mt = mt; c.xrt = xrt; c.xrt_cl = xrt; c.m_chord = m_chord; c.tchord = tchord;


%% write file
directory = job.rjm_directory;
if write_file == 1	
%% Write blade file
fid = fopen([directory 'blade.mises'],'w');
fprintf(fid,'%s\n','MISES_Section');
fprintf(fid,'%6.5f %6.5f %3.2f %3.2f %7.6f\n',...
    [tand(c.chi_le) tand(c.chi_te) m_chord m_chord 2*pi / c.N]);
% mt(:,2) = - mt(:,2);
for i = 1:size(mt,1)
    fprintf(fid,'%10.12f %10.12f\n',mt(i,:));
end
fclose(fid);
end

end