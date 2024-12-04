function Blade=ReadBlade(varargin);
% Function to read Mises blade.xxx files.
% --
% Initialise variables
switch nargin; % Deal with inputs
   case {0}; F.ext='mises'; F.path=pwd; % Assume current directory and file extension 'mises'
   case {1}; F.ext=varargin{1}; F.path=pwd;
   case {2}; F.ext=varargin{1}; F.path=varargin{2};
end
clear varargin
%fprintf(' Reading blade file...\n');

% Generate name and open blade file for reading
F.filename=['blade.',F.ext]; F.file=fullfile(F.path,F.filename); F.fid=fopen(F.file,'r');
if F.fid < 0; % Return error code if there is a problem
   warning(['Problem opening "',F.filename,'". Returning error code.']);
   Blade=-1;
else % Read in the file
   % Read title
   Blade.title = fgetl(F.fid);
   % Read header line
   line = fgetl(F.fid);
   store = sscanf(line,'%f');
   Blade.sinl = store(1);
   Blade.sout = store(2);
   Blade.chinl = store(3);
   Blade.chout = store(4);
   Blade.pitch = store(5);
   % Read blade geometry
   ib = 0; %blade number
   nextBlade = 1;
   while nextBlade > 0
      ib = ib + 1;
      ii = 0;
      nextpt = 1;
      while nextpt > 0
         ii = ii + 1;
         line = fgetl(F.fid);
         % HOT: Changed this section as original led to the possibility of
         % the last point being missed out. This change ensures that this
         % never happens.
         %if feof(F.fid) | isempty(line)
         if isempty(line)
            nextpt = 0;
            nextBlade = 0;
         else
            store = sscanf(line,'%f');
            if store(1) == 999
               nextpt = 0;
            else
               Blade.x{ib}(ii) = store(1);
               Blade.y{ib}(ii) = store(2);
            end
         end
         if feof(F.fid)
             nextpt = 0;
             nextBlade = 0;
         end
         % HOT: End of changes
      end
   end
   Blade.nbld = ib;
   % Order struct
   Blade=orderfields(Blade);
   % Close blade file
   t.isclosed=fclose(F.fid);
   if t.isclosed < 0; warning(['Problem closing "',F.filename,'"']); end
   clear t
end

return