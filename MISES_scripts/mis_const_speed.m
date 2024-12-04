function [job]=mis_const_speed(directory)
% function to run along a constant speed line (create new ises file and run
% polar), after running mis_run_section

is_turb = 1;

directory = strrep(directory,'TURBOSTREAM','MISES');

% Read in flow file
if exist([directory 'polarx.mises'],'file') ~= 0
    [Polarx, Ises] = mis_read_polarx('mises',directory);
else
    disp('File Not Found')
    p = [];
    return
end

[~,j] = min(abs(Polarx.binl-Ises.binl));

load([directory 'section.mat']);

f.M = Polarx.minl(j);
f.M2 = Polarx.mout(j);
f.Alpha = Polarx.binl(j);
f.Alpha2 = Polarx.bout(j);
f.p1pt = Polarx.p1pt(j);
f.p2pt = Polarx.p2pt(j);
% calculate Vta01
ga=1.4;
f.Vtao1 = f.M/(sqrt(1+(ga-1)/2*f.M^2))*tand(abs(f.Alpha))/(sqrt(1+tand(abs(f.Alpha))^2));

% Maximum MISES length
m_chord = max(c.mt(:,1));

% create new directory
%% make directory
% 	directory_smooth = [directory(1:end-1) num2str(n) '/'];
directory_conspeed = [directory(1:end-1) 'f/'];

job.dir_conspeed = directory_conspeed;

% Make directory if required
if exist(directory_conspeed,'dir') == 0
	mkdir(directory_conspeed);
end

[~,~]=system(['cp -r ' directory  '* ' directory_conspeed]);
[~,~]=system(['rm ' directory_conspeed 'polarx.mises'],'-echo');


% save new section.mat
save([directory_conspeed 'section.mat'],'c','f');

% create new ises file
%% change ises back to normal smoothing
% Write ises file
fid = fopen([directory_conspeed 'ises.mises'],'w');
% fprintf(fid,'%s\n','1 2 5  !Global variables');
% fprintf(fid,'%s\n','1 3 4  !Global constraints');
% fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
% fprintf(fid,'%s\n',' 3 4 6 9 17 !Global constraints');
fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
fprintf(fid,'%s\n',' 3 4 6 9 17 !Global constraints');


fprintf(fid,'%6.5f %6.5f %6.5f %6.5f %6.5f ',[f.M f.p1pt tand(f.Alpha) -0.66*m_chord f.Vtao1]);
% fprintf(fid,'%6.5f %6.5f %6.5f %6.5f %6.5f ',[0.75 f.Pin/f.Poin tand(f.Alpha) -0.66*m_chord f.Vtao1]);
fprintf(fid,'%s\n','|Minl p1/p01 Sinl Xinl v1/ao1');
% fprintf(fid,'0.00000 0.00000 0.00000 %6.5f |Mout p2/p01 Sout Xout\n',1.33*m_chord);
fprintf(fid,'%6.5f %6.5f %6.5f %6.5f  ', [f.M2 f.p2pt tand(f.Alpha2) 1.33*m_chord]);
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

cd(directory_conspeed)
[~,~]=system(['rm ' directory_conspeed 'polarx.mises'],'-echo');
[~,~]=system(['rm ' directory_conspeed 'idat.mises*']);
[n_type]=mis_run_inner(directory_conspeed,17);

% run polar

% Write ises file
fid = fopen([directory_conspeed 'ises.mises'],'w');
% fprintf(fid,'%s\n','1 2 5  !Global variables');
% fprintf(fid,'%s\n','1 3 4  !Global constraints');
% fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
% fprintf(fid,'%s\n',' 3 4 6 9 17 !Global constraints');
fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
fprintf(fid,'%s\n',' 3 4 6 9 17 !Global constraints');


fprintf(fid,'%6.5f %6.5f %6.5f %6.5f %6.5f ',[f.M f.p1pt tand(f.Alpha) -0.66*m_chord f.Vtao1]);
% fprintf(fid,'%6.5f %6.5f %6.5f %6.5f %6.5f ',[0.75 f.Pin/f.Poin tand(f.Alpha) -0.66*m_chord f.Vtao1]);
fprintf(fid,'%s\n','|Minl p1/p01 Sinl Xinl v1/ao1');
% fprintf(fid,'0.00000 0.00000 0.00000 %6.5f |Mout p2/p01 Sout Xout\n',1.33*m_chord);
fprintf(fid,'%6.5f %6.5f %6.5f %6.5f  ', [f.M2 f.p2pt tand(f.Alpha2) 1.33*m_chord]);
fprintf(fid,'%s\n','|Mout p2/p01 Sout Xout');
fprintf(fid,'%s\n','0.00000     0.00000     1.39433     0.00000  | MFR  HWRATin GAMin DIFACTOR');
fprintf(fid,'%i ',round(f.Re));
fprintf(fid,'%s\n', '-4.00000 0.00000  0.26133 |Re -Turb');
if is_turb == 1
    fprintf(fid,'%s\n','0.001 0.01  |Xtr1 Xtr2');
else
    fprintf(fid,'%s\n','1.10 1.10  |Xtr1 Xtr2');
end
fprintf(fid,'%s\n','4     0.95000  -1.00000  0.00000   | ISMOM  MCRIT  MUCON  PLOSSIN');
fprintf(fid,'%s\n','0.00000  0.00000   | BVR1in  BVR2in');
fclose(fid);

[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMises_Ho-on2.py"','-echo');	

%% check if you need to run polar again to get a fuller loss loop
% Get all idat filenames
A = [dir([directory_conspeed 'idat.mises_0*'])];
F = cell(length(A),1); for n = 1:length(A); F{n} = A(n).name; end;

% check if you need to run again
% Get all idat filenames
A2 = [dir([directory_conspeed 'idat.mises_5*'])];
F2 = cell(length(A2),1); for n = 1:length(A2); F2{n} = A2(n).name; end;

F
F2
if length(F) > 4 && length(F2) > 4
else
	
directory_conspeed = [directory(1:end-1) 'f2/'];
% Make directory if required
if exist(directory_conspeed,'dir') == 0
	mkdir(directory_conspeed);
end

[~,~]=system(['cp -r ' job.dir_conspeed  '* ' directory_conspeed],'-echo');
[~,~]=system(['rm ' directory_conspeed '*polar*'],'-echo');
[~,~]=system(['rm ' directory_conspeed '*idat*'],'-echo');
	
%% change ises back to normal smoothing
% Write ises file
fid = fopen([directory_conspeed 'ises.mises'],'w');
% fprintf(fid,'%s\n','1 2 5  !Global variables');
% fprintf(fid,'%s\n','1 3 4  !Global constraints');
% fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
% fprintf(fid,'%s\n',' 3 4 6 9 17 !Global constraints');
fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
fprintf(fid,'%s\n',' 3 4 6 9 17 !Global constraints');


fprintf(fid,'%6.5f %6.5f %6.5f %6.5f %6.5f ',[f.M f.p1pt tand(f.Alpha) -0.66*m_chord f.Vtao1]);
% fprintf(fid,'%6.5f %6.5f %6.5f %6.5f %6.5f ',[0.75 f.Pin/f.Poin tand(f.Alpha) -0.66*m_chord f.Vtao1]);
fprintf(fid,'%s\n','|Minl p1/p01 Sinl Xinl v1/ao1');
% fprintf(fid,'0.00000 0.00000 0.00000 %6.5f |Mout p2/p01 Sout Xout\n',1.33*m_chord);
fprintf(fid,'%6.5f %6.5f %6.5f %6.5f  ', [f.M2 f.p2pt tand(f.Alpha2) 1.33*m_chord]);
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

cd(directory_conspeed)
[~,~]=system(['rm ' directory_conspeed 'polarx.mises']);
[~,~]=system(['rm ' directory_conspeed 'idat.mises*']);
[n_type]=mis_run_inner(directory_conspeed,17,n_type);

% run polar

% Write ises file
fid = fopen([directory_conspeed 'ises.mises'],'w');
% fprintf(fid,'%s\n','1 2 5  !Global variables');
% fprintf(fid,'%s\n','1 3 4  !Global constraints');
% fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
% fprintf(fid,'%s\n',' 3 4 6 9 17 !Global constraints');
fprintf(fid,'%s\n',' 1 2 5 6 15 !Global variables');
fprintf(fid,'%s\n',' 3 4 6 9 17 !Global constraints');


fprintf(fid,'%6.5f %6.5f %6.5f %6.5f %6.5f ',[f.M f.p1pt tand(f.Alpha) -0.66*m_chord f.Vtao1]);
% fprintf(fid,'%6.5f %6.5f %6.5f %6.5f %6.5f ',[0.75 f.Pin/f.Poin tand(f.Alpha) -0.66*m_chord f.Vtao1]);
fprintf(fid,'%s\n','|Minl p1/p01 Sinl Xinl v1/ao1');
% fprintf(fid,'0.00000 0.00000 0.00000 %6.5f |Mout p2/p01 Sout Xout\n',1.33*m_chord);
fprintf(fid,'%6.5f %6.5f %6.5f %6.5f  ', [f.M2 f.p2pt tand(f.Alpha2) 1.33*m_chord]);
fprintf(fid,'%s\n','|Mout p2/p01 Sout Xout');
fprintf(fid,'%s\n','0.00000     0.00000     1.39433     0.00000  | MFR  HWRATin GAMin DIFACTOR');
fprintf(fid,'%i ',round(f.Re));
fprintf(fid,'%s\n', '-4.00000 0.00000  0.26133 |Re -Turb');
if is_turb == 1
    fprintf(fid,'%s\n','0.001 0.01  |Xtr1 Xtr2');
else
    fprintf(fid,'%s\n','1.10 1.10  |Xtr1 Xtr2');
end
fprintf(fid,'%s\n','4     0.95000  -1.00000  0.00000   | ISMOM  MCRIT  MUCON  PLOSSIN');
fprintf(fid,'%s\n','0.00000  0.00000   | BVR1in  BVR2in');
fclose(fid);

[~,~] = system('bash --login -c "python ~/bin/CalculateBladeMises_Ho-on2.py"','-echo');	

[~,~]=system(['cp -r ' directory_conspeed  'idat.mises_* ' job.dir_conspeed]);
end

end