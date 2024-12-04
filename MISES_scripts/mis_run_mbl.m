function p = mis_run_mbl(directory,c,f,h,C,plot_stuff)
% MIS_RUN_MBL  Run mises boundary layer solver on a specified pressure distribution
%
%   p = MIS_RUN_MBL(directory,c,f,h,C,plot_stuff)
%
%   directory - string of folder to run in
%   c - data structure containing coordinates and static pressure distribution
%   f - data structure of inlet boundary conditions
%   h - optional figure handle
%   C - optional RGB colour array
%   plot_stuff - 0 or 1 for showing working
%   p - output data structure

if exist('C','var') == 0 || isempty(C) == 1
    C = [0 0 0];
end
if exist('plot_stuff','var') == 0
    plot_stuff = 1;
end
if exist('h','var') == 0 || isempty(h) == 1 && plot_stuff == 1
    h = figure();
end

% Default to zero contraction
if isfield(c,'cr') == 0
    cr = ones(size(c.x));
end

% Find stagnation point
[~,i_stag] = max(c.P);

% Find 98% chord position
t = parameterise_section(c.x,c.rt);
[~,i_te] = min(abs(t.s_2_raw - 0.98));
i_te = i_te+i_stag;

% Write streamwise coordinate
s = zeros(size(c.x));
s(t.i_1) = t.s_1_raw;
s(t.i_2) = t.s_2_raw;
s = s(i_stag:i_te); 
s = [0 ; cumsum(abs(diff(s)))];
s = (s - min(s)) / (max(s) - min(s));

% Calculate isentropic edge velocity
Ue = abs(((f.Po - c.P(i_stag:i_te)) / (0.5 * f.ro)).^0.5);

% Smooth distributions of points and edge velocity
knots = [0 0 hyperbolic_bunch(30,0.004,0.001) 1 1];
s_interp = hyperbolic_bunch(200,1e-3,1e-3).';

spline_temp = spap2(knots,3,s,Ue);
Ue = fnval(spline_temp,s_interp);

spline_temp = spap2(knots,3,s,c.r(i_stag:i_te));
r = fnval(spline_temp,s_interp);

spline_temp = spap2(knots,3,s,c.cr(i_stag:i_te));
cr = fnval(spline_temp,s_interp);

s = s_interp;

% Noramlise edge velocity
Ue = Ue/((f.T*1.4*287)^0.5);

% Plot input
% figure(h)
% plot(s,Ue)

% Delete previous files
A = dir([directory '*']);
for a = length(A)
    delete([directory A(a).name]);
end

% Write input file
% [M Re rot]
% Hwrat
% [Xtrip Ncrit]
% [ x Ue r cr]
fid = fopen([directory 'input.mblrun.mises'],'w');

fprintf(fid,'%5.4f\t%10.2f\t%5.4f\t| M Re Rot\n',f.M,f.Re,f.rpm);
fprintf(fid,'0.0000\t| Hwrat\n');
fprintf(fid,'%6.4f\t%5.4f\t| Xtrip Ncrit\n',0.02,9);
for n = 1:length(s)
    fprintf(fid,'%7.5f\t%7.5f\t%7.5f\t%7.5f\n',s(n),Ue(n),r(n)/t.chord,cr(n));
end
fclose(fid);

% Run boundary layer solver
cd(directory)
setenv('GFORTRAN_STDIN_UNIT', '5') 
setenv('GFORTRAN_STDOUT_UNIT', '6') 
setenv('GFORTRAN_STDERR_UNIT', '0')
[~,~] = system('bash --login -c "python ~/bin/CalculateBoundaryLayer.py"');
setenv('GFORTRAN_STDIN_UNIT', '-1') 
setenv('GFORTRAN_STDOUT_UNIT', '-1') 
setenv('GFORTRAN_STDERR_UNIT', '-1')
cd('~/Documents/MATLAB/')

% Read output file and plot boundary layer parameters
A = dlmread([directory 'mblrun.mises'],'',5,0);

if plot_stuff == 1
    figure(h);
    subplot(2,2,1); hold on; grid on; box on;
    plot(A(:,1),A(:,2),'-','Color',C)
    xlabel('Chord'); ylabel('Edge Velocity');

    subplot(2,2,2); hold on; grid on; box on;
    plot(s,cr,'-','Color',C)
    xlabel('Chord'); ylabel('Contraction Ratio');

    subplot(2,2,3); hold on; grid on; box on;
    plot(A(:,1),A(:,6),'-','Color',C)
    xlabel('Chord'); ylabel('Momentum Thickness');

    subplot(2,2,4); hold on; grid on; box on;
    plot(A(:,1),A(:,7),'-','Color',C)
    xlabel('Chord'); ylabel('Shape Factor');
end

p.s = A(:,1);
p.Ue = A(:,2);
p.theta = A(:,6);
p.H = A(:,7);
p.cr = cr;

end