function mis_write_input(directory,c,f,is_turb)
% MIS_WRITE_INPUT  Write input files for a given section geometry 
%
%   MIS_WRITE_INPUT(directory,c,f,is_turb)
%
%   directory - string of directory to run in
%   c - input section data structure or coordinates
%   f - input boundary conditions
%   is_turb - optional transition inputs
%
%   c is either data structure or cell array:
%       if data structure then it must be a full definition to run in
%       bl_construct_section
%       if cell array then cell 1 is xrrt coords, cell 2 is blade count


% Make directory if required
if exist(directory,'dir') == 0
    mkdir(directory);
end

% Delete old output files from directory
if isempty(dir([directory 'polar*'])) == 0; delete([directory 'polar*']); end;
if isempty(dir([directory 'output*'])) == 0; delete([directory 'output*']); end;
if isempty(dir([directory 'idat*'])) == 0; delete([directory 'idat*']); end;
if isempty(dir([directory 'spec*'])) == 0; delete([directory 'spec*']); end;


%% Prepare coordinates for MISES

% Determine whether to run on a section definition or directly on
% coordinates
if isstruct(c) == 1
    
    % Shift centroid to origin
    c.xy_cen = [0 0];
    
    % Construct blade from spline parameters
    [xrt_cl,xrt_cam] = bl_construct_section(c,0,0);
    xrt = bl_construct_section(c,0,1);

    % Record radii at leading and trailing edge
    r_le = c.r_le; r_te = c.r_te;
    
    % Record camber line coordinates
    c.xrt_cam = xrt_cam;
    
elseif iscell(c) == 1
    
    % Run directly from specified geometry of closed section
    xrt_cl = c{1}(:,[1 3]); N = c{2};
    
    % Parameterise to find trailing edge circle
    t = bl_parameterise_section(xrt_cl(:,1:2)); 
    
    % Cut out open trailing edge section
    xrt = [xrt_cl(max(t.i_circ):end,:) ; xrt_cl(2:min(t.i_circ),:)];
    
    % Record radii at leading and trailing edge
    r_le = c{1}(t.i_1(1),2); r_te = c{1}(t.i_1(end),2); 
    
    % Shift to origin
    xrt(:,1) = xrt(:,1) - min(xrt(:,1));
    
    % Save data structure of section parameters
%     r_le = c{1}(t.i_le,2); 
    clear c; c.N = N; c.tchord = t.tchord; c.chi_le = t.chi_le; c.chi_te = t.chi_te;
    c.xrt_cam = t.xy_cam; c.r_le = r_le; c.r_te = r_te;
    c.thick_max = t.thick_max; c.thick_te = t.thick_te; c.dcam_le = t.dcam_le; c.dcam_te = t.dcam_te;
    c.s_thick_max = t.s_thick_max;
end

% Calculate radial coordinates
r = interp1(c.xrt_cam([1 end],1),[c.r_le c.r_te],xrt(:,1),'linear','extrap');
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
if f.Alpha < 0; f.Alpha = -f.Alpha; end;

% Calculate MISES coordinates
mt = [x ./ r rt ./ r]; 

% Maximum MISES length
m_chord = max(mt(:,1));

% Record leading edge point for plotting later
c.mt_le = c.xrt_cam(1,:) / c.r_le;

% Parameterise the camber line
t = bl_parameterise_camber(c.xrt_cam); c.s_cl = t.s_cl; c.cam = t.cam;

% Interpolate the leading edge centre for plotting later
xrt_le_cen = [interp1(c.s_cl,c.xrt_cam(:,1),0.03,'pchip') interp1(c.s_cl,c.xrt_cam(:,2),0.03,'pchip')];
c.mt_le_cen = xrt_le_cen / interp1(c.xrt_cam([1 end],1),[c.r_le c.r_te],xrt_le_cen(1));

% Save remaining coordinates to section data structure
c.mt = mt; c.xrt = xrt; c.xrt_cl = xrt_cl;


%% Write input files

% Write section file
save([directory 'section.mat'],'c','f');

% Write blade file
fid = fopen(fullfile(directory, 'blade.mises'),'w');
fprintf(fid,'%s\n','MISES_Section');
fprintf(fid,'%6.5f %6.5f %3.2f %3.2f %7.6f\n',...
    [tand(c.chi_le) tand(c.chi_te) m_chord m_chord 2*pi / c.N]);
% mt(:,2) = - mt(:,2);
for i = 1:size(mt,1)
    fprintf(fid,'%10.12f %10.12f\n',mt(i,:));
end
fclose(fid);

% Write stream file
fid = fopen(fullfile(directory, 'stream.mises'),'w');
fprintf(fid,'%i\n',0);
fprintf(fid,'%6.5f %6.5f %6.5f\n',[min(mt(:,1))-0.5 r_le/c.tchord 1.00000]);
fprintf(fid,'%6.5f %6.5f %6.5f\n',[0.0 r_le/c.tchord 1.00000]);
fprintf(fid,'%6.5f %6.5f %6.5f\n',[0.0 r_le/c.tchord 1.00000]);
fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1)) r_te/c.tchord f.AVDR]);
fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1)) r_te/c.tchord f.AVDR]);
fprintf(fid,'%6.5f %6.5f %6.5f\n',[max(mt(:,1))+0.5 r_te/c.tchord f.AVDR]);
fclose(fid);

% Write ises file
fid = fopen(fullfile(directory, 'ises.mises'),'w');
fprintf(fid,'%s\n','1 2 5  !Global variables');
fprintf(fid,'%s\n','1 3 4  !Global constraints');
fprintf(fid,'%6.5f %6.5f %6.5f %6.5f  ',[f.M 0.0 tand(f.Alpha) -0.66*m_chord]);
fprintf(fid,'%s\n','|Minl p1/p01 Sinl Xinl');
fprintf(fid,'0.00000 0.00000 0.00000 %6.5f |Mout p2/p01 Sout Xout\n',1.33*m_chord);
fprintf(fid,'%s\n','0.00000     0.00000     1.39433     0.00000  | MFR  HWRATin GAMin DIFACTOR');
fprintf(fid,'%i ',round(f.Re));
fprintf(fid,'%s\n', '-2.00000 0.00000  0.26133 |Re -Turb');
if is_turb == 1
    fprintf(fid,'%s\n','0.02 0.02  |Xtr1 Xtr2');
else
    fprintf(fid,'%s\n','1.10 1.10  |Xtr1 Xtr2');
end
fprintf(fid,'%s\n','4     0.95000  1.00000  0.00000   | ISMOM  MCRIT  MUCON  PLOSSIN');
fprintf(fid,'%s\n','0.00000  0.00000   | BVR1in  BVR2in');
fprintf(fid,'%s\n','5.6  1.0  1.0  | SCC SCP SCD');
fprintf(fid,'%s\n','-0.071652  1.0  | XSHOCK FCTSHK');

fclose(fid);

% Write gridpar file
fid = fopen(fullfile(directory, 'gridpar.mises'),'w');
fprintf(fid,'%s\n','T T');
fprintf(fid,'%s\n','30   30');
fprintf(fid,'%s\n','24');
fprintf(fid,'%s\n','1.500000');
fprintf(fid,'%s\n','0.800000');
fprintf(fid,'%s\n','80    0.10000    0.900000    1.000000');
fprintf(fid,'%s\n','1.000000    1.000000    0.000000    1.000000    1.000000    0.000000');
fclose(fid);

% Run iset to initialise the grid and solution
mis_run(directory,'init');


end

