function Ises = mis_read_ises(varargin)
% Function to read Mises ises.xxx files.
% --
% Initialise variables
switch nargin; % Deal with inputs
   case {0}; F.ext='mises'; F.path=pwd; % Assume current directory and file extension 'mises'
   case {1}; F.ext=varargin{1}; F.path=pwd;
   case {2}; F.ext=varargin{1}; F.path=varargin{2};
end
clear varargin
%fprintf(' Reading ises file...\n');

% Generate name and open ises file for reading
F.filename=['ises.',F.ext]; F.file=fullfile(F.path,F.filename); F.fid=fopen(F.file,'r');
if F.fid < 0; % Return error code if there is a problem
   warning(['Problem opening "',F.filename,'". Returning error code.']);
   Ises=-1;
else % Read in the file
   % gvar
   line = fgetl(F.fid);
   Ises.gvar = sscanf(line,'%i');
   % gconccc
   line = fgetl(F.fid);
   Ises.gcon = sscanf(line,'%i');
   % inlet bcs
   line = fgetl(F.fid);
   store = sscanf(line,'%f');
   Ises.minl = store(1);
   Ises.p1pt = store(2);
   Ises.sinl = store(3);
   Ises.xinl = store(4);
   if length(store) == 5
      Ises.v1at = store(5);
   end
   % outlet bcs
   line = fgetl(F.fid);
   store = sscanf(line,'%f');
   Ises.mout = store(1);
   Ises.p2pt = store(2);
   Ises.sout = store(3);
   Ises.xout = store(4);
   if length(store) == 5
      Ises.v2at = store(5);
   end
   % Splitter mass frac and non-adiabatic wall
   line = fgetl(F.fid);
   store = sscanf(line,'%f');
   Ises.mfr = store(1);
   Ises.hwrat = store(2);
   Ises.gamin = store(3);
   % Reynolds number and turbulence levels
   line = fgetl(F.fid);
   store = sscanf(line,'%f');
   Ises.reyn = store(1);
   Ises.ncrit = store(2);
   % Forced transition points
   line = fgetl(F.fid);
   store = sscanf(line,'%f');
   Ises.strp = [store(1), store(2)];
   % Algorithm controls
   line = fgetl(F.fid);
   store = sscanf(line,'%f');
   Ises.ismom = store(1);
   Ises.mcrit = store(2);
   Ises.mucon = store(3);
   % Stream-tube thickness mode amplitudes
   if ~feof(F.fid)
      line = fgetl(F.fid);
      store = sscanf(line,'%f');
      Ises.bvr=[store(1),store(2)];
   end
   if ~feof(F.fid)
         % Shear lag consts
       line = fgetl(F.fid);
       store = sscanf(line,'%f');
       Ises.Klag = store(1);
       Ises.Kp = store(2);
       Ises.Kd= store(3);
       % Ctau source consts
       line = fgetl(F.fid);
       store = sscanf(line,'%f');
       Ises.Xshock = store(1);
       Ises.Fctshk = store(2);
   end
   % Geometry movement, scaling and rotation mode amplitudes
   if ~feof(F.fid)
      line = fgetl(F.fid);
      store = sscanf(line,'%f');
      Ises.movx = store(1);
      Ises.movy = store(2);
      Ises.scal = store(3);
      Ises.rota = store(4);
   end
   % Geometry shape mode amplitudes
   if ~feof(F.fid)
      line = fgetl(F.fid);
      store = sscanf(line,'%f');
      Ises.kmod = store(1);
      Ises.gmod = store(2);
   end
   % Deal with turbulence intensity
   if Ises.ncrit>=0 % Just use value for e^n transition model
      Ises.ncrit=Ises.ncrit;
      t.TuPrime=(100.*(exp((8.43+Ises.ncrit)./(-2.4))));
      Ises.turb=((2.7./2).*(log((1+(t.TuPrime./2.7))./(1-(t.TuPrime./2.7)))));
   elseif Ises.ncrit<0 % Calculate the value of ncrit used by mises (see manuals)
      Ises.turb=-Ises.ncrit;
      t.TuPrime=(2.7.*tanh((Ises.turb)./2.7));
      Ises.ncrit=(-8.43-(2.4.*log(t.TuPrime./100)));
   end
   % Calculate flow angles for later use
   Ises.binl=((180/pi)*atan(Ises.sinl));
   Ises.bout=((180/pi)*atan(Ises.sout));
   % Order struct
   Ises=orderfields(Ises);
   % Close idat.* file
   t.isclosed=fclose(F.fid);
   if t.isclosed < 0; warning(['Problem closing "',F.filename,'"']); end
   clear t
end

return