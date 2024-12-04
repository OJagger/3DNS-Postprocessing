function Stream=ReadStream(varargin)
% Function to read Mises stream.xxx files.
% --
% Initialise variables
switch nargin; % Deal with inputs
   case {0}; F.ext='mises'; F.path=pwd; % Assume current directory and file extension 'mises'
   case {1}; F.ext=varargin{1}; F.path=pwd;
   case {2}; F.ext=varargin{1}; F.path=varargin{2};
end
clear varargin
%fprintf(' Reading stream file...\n');

% Generate name and open stream file for reading
F.filename=['stream.',F.ext]; F.file=fullfile(F.path,F.filename); F.fid=fopen(F.file,'r');
if F.fid < 0; % Return error code if there is a problem
   warning(['Problem opening "',F.filename,'". Returning error code.']);
   Stream=-1;
else % Read in the file
   % Read in rotation speed
   line=fgetl(F.fid);
   Stream.rotrel=sscanf(line,'%f',1);
   % Read streamtube coordinates
   ii=1; carryon=1;
   while carryon > 0
      line = fgetl(F.fid);
      if line==-1
         carryon = 0;
      else
         store = sscanf(line,'%f');
         Stream.x(ii,1) = store(1);
         Stream.r(ii,1) = store(2);
         Stream.b(ii,:) = store(3:end);
      end
      ii = ii + 1;
   end
   % Delete two duplicated points at le and te
%    for ii=2:length(Stream.x)-2;
%        ii
%       if Stream.x(ii,1)==Stream.x(ii-1,1);
%          Stream.x(ii,:)=[]; Stream.r(ii,:)=[]; Stream.b(ii,:)=[];
%       end
%    end
    % delete duplicate points
    [Stream.x,i_uniq] = unique(Stream.x);
    Stream.r = Stream.r(i_uniq);
    Stream.b = Stream.b(i_uniq);


   % Order struct
   Stream=orderfields(Stream);
   % Close stream output file
   t.isclosed=fclose(F.fid);
   if t.isclosed < 0; warning(['Problem closing "',F.filename,'"']); end
   clear t
end

return