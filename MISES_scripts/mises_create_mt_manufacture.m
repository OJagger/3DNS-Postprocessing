% create mt from job.mat file 
% option to write out blade file

function [mt,c] = mises_create_mt(job,write_file,c)

if exist('write_file','var')==0
	write_file = 1; % write blade.mises file by default
end

directory = strrep(job.rjm_directory,'TURBOSTREAM','MISES');

[xrt,sl_cx,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
% Distance through surfaces
s_1 = [0 ; cumsum(sum(diff(mca.xrt_ss,1,1).^2,2).^0.5)]; s_1 = s_1 / max(s_1);
s_2 = [0 ; cumsum(sum(diff(mca.xrt_ps,1,1).^2,2).^0.5)]; s_2 = s_2 / max(s_2);
% remove trailing edge
i_circ = bl_find_circ(mca.xrt_ss,mca.xrt_ps,s_1,s_2);
% Cut out open trailing edge section
xrt = [flip(mca.xrt_ss(1:i_circ,:)) ; mca.xrt_ps(2:i_circ(2),:)];

r = interp1(mca.xrt_cam([1 end],1),[job.r_mid job.r_mid],xrt(:,1),'linear','extrap');
xrrt = [xrt(:,1) r xrt(:,2)];

% Downsample section
s = dist_2d(xrrt,1);
s_interp = interp1(linspace(0,1,length(s))',s,linspace(0,1,250)');
% i = round(linspace(1,size(xrrt,1),250));
% x = xrrt(i,1); r = xrrt(i,2); rt = xrrt(i,3);
x = interp1(s,xrrt(:,1),s_interp,'pchip');
r = interp1(s,xrrt(:,2),s_interp,'pchip');
rt = interp1(s,xrrt(:,3),s_interp,'pchip');



% Sample equal spacing of points
% s = dist_2d(xrrt,1);
% x = interp1(s,xrrt(:,1),linspace(0,1,250)');
% r = interp1(s,xrrt(:,2),linspace(0,1,250)');
% rt = interp1(s,xrrt(:,3),linspace(0,1,250)');

% Check blade nose always points down and inlet flow angle is positive
[~,i_min] = min(x);
if rt(i_min) > rt(1); rt = -rt; end; 


% Calculate MISES coordinates
mt = [x ./ r rt ./ r]; 

% Maximum MISES length
m_chord = max(mt(:,1));

% Record leading edge point for plotting later
% c.mt_le = c.xrt_cam(1,:) / c.r_le;
c.mt_le = mca.xrt_cam(1,:) / c.r_le;

%% need to re-do this
% Parameterise the camber line
t = bl_parameterise_camber(c.xrt_cam); c.s_cl = t.s_cl; c.cam = t.cam;

% Interpolate the leading edge centre for plotting later
xrt_le_cen = [interp1(c.s_cl,c.xrt_cam(:,1),0.03,'pchip') interp1(c.s_cl,c.xrt_cam(:,2),0.03,'pchip')];
c.mt_le_cen = xrt_le_cen / interp1(c.xrt_cam([1 end],1),[c.r_le c.r_te],xrt_le_cen(1));
%%

% Save remaining coordinates to section data structure
% c.mt = mt; c.xrt = xrt; c.xrt_cl = xrt_cl;
c.mt = mt; c.xrt = xrt; c.xrt_cl = xrt;

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