function [Polar]=mis_read_polar(varargin)
% Function to read Mises polar.xxx files.
% --
% Initialise variables
switch nargin; % Deal with inputs
   case {0}; F.ext='mises'; F.path=pwd; % Assume current directory and file extension 'mises'
   case {1}; F.ext=varargin{1}; F.path=pwd;
   case {2}; F.ext=varargin{1}; F.path=varargin{2};
end

% keep using job.rjm_directory
F.path = strrep(F.path,'TURBOSTREAM','MISES');
cd(F.path); clear varargin


fprintf(' Reading polar file...\n');

% Generates name and opens polar file for reading
F.filename=['polar.',F.ext]; F.file=fullfile(F.path,F.filename); F.fid=fopen(F.file,'r');

if F.fid < 0; % Return error code if there is a problem
   warning(['Problem opening "',F.filename,'". Returning error code.']);
   Polar=-1;
else % Read in the file
   % Skip 1st 7 lines and read title
   t.dd=fgetl(F.fid);
   t.dd=fgetl(F.fid);
   t.dd=fgetl(F.fid);
   t.title=fgetl(F.fid);
   Polar.title=t.title(24:56);
   Polar.title=deblank(Polar.title);
   t.dd=fgetl(F.fid);
   t.dd=fgetl(F.fid);
   t.dd=fgetl(F.fid);
   clear t;
   % Read polar data
   carryon=1; ii=0;
   while carryon > 0
      ii=ii+1;
      t.line=fgetl(F.fid);
      if ischar(t.line)
         t.line=sscanf(t.line,'%f');
         if ~any(isnan(t.line))
            ainl(ii,1)=atan(t.line(1)).*180./pi;
            aout(ii,1)=atan(t.line(2)).*180./pi;
            minl(ii,1)=t.line(3);
            mout(ii,1)=t.line(4);
            pinlp0(ii,1)=t.line(5);
            poutp0(ii,1)=t.line(6);
            re(ii,1)=t.line(7).*1e6;
            tu(ii,1)=t.line(8)./100;
            yp(ii,1)=t.line(9);
            ypv(ii,1)=t.line(10);
            xtr1(ii,1)=t.line(11);
            xtr2(ii,1)=t.line(12);
            drvt(ii,1)=t.line(13);
         else
            ii=ii-1;
         end
      else
         carryon = 0;
      end
   end
   % Convergence failure flag
   Polar.nsoln = ii-1;
   if ii <= 1
       Polar = -1;
       return
   end
   if isstruct(Polar) % Reorder polar into ascending order if it's a struct
      t.dat=[ainl aout minl mout pinlp0 poutp0 re tu yp ypv xtr1 xtr2 drvt];
      t.datnew=sortrows(t.dat,1);
      Polar.ainl=t.datnew(:,1);
      Polar.aout=t.datnew(:,2);
      Polar.minl=t.datnew(:,3);
      Polar.mout=t.datnew(:,4);
      Polar.pinlp0=t.datnew(:,5);
      Polar.poutp0=t.datnew(:,6);
      Polar.re=t.datnew(:,7);
      Polar.tu=t.datnew(:,8);
      Polar.yp=t.datnew(:,9);
      Polar.ypv=t.datnew(:,10);
      Polar.xtr1=t.datnew(:,11);
      Polar.xtr2=t.datnew(:,12);
      Polar.drvt=t.datnew(:,13);
   end
   % Order struct
   Polar=orderfields(Polar);
   % Close idat.* file
   t.isclosed=fclose(F.fid);
   if t.isclosed < 0; warning(['Problem closing "',F.filename,'"']); end
   clear t
end

return