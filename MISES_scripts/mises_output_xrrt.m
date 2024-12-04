% create xrrt (c structure in mises_run_section) from job.mat file % dl467
function [c]=mises_output_xrrt(job)

% section radius
if isfield(job,'r_le')==0
   r_le = job.r_mid; 
else
    r_le =job.r_mid;
end
if isfield(job,'r_te')==0
   r_te = job.r_mid; 
else
    r_te =job.r_mid;
end
c.r_le = r_e;
c.r_te= r_te;

if isfield(job,'N') == 1
    c.N = job.N;
end

% calculate pitch
pitch = 2*pi*0.5*(c.r_le+c.r_te)./c.N;


ote = 0; dev =[]; plot_stuff = 0;
cam_spline = 1;
% construct section with trailing edge
[xrt,xrt_cam,xrt_ps,xrt_ss]=bl_construct_section(job,plot_stuff,ote,dev,cam_spline);

[c,h] = bl_parameterise_section(xrt,0,0,pitch);

% copy new fields into c
c.r_le = r_le;
c.r_te= r_te;
c.N = job.N;

end