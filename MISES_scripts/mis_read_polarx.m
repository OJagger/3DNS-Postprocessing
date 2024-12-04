function [Polarx, Ises, Blade, Stream, raw] = mis_read_polarx(varargin)
% Function to read Mises polarx.xxx files.
% --
% Initialise variables
switch nargin; % Deal with inputs
   case {0}; F.ext='mises'; F.path=pwd; % Assume current directory and file extension 'mises'
   case {1}; F.ext=varargin{1}; F.path=pwd;
   case {2}; F.ext=varargin{1}; F.path=varargin{2};
end
%cd(F.path)
F.path = strrep(F.path,'TURBOSTREAM','MISES');
clear varargin
F.nsmax=100; % Set maximum number of solutions that can be read in
%fprintf(' Reading polarx file...\n');

% Read Ises.*, Blade.* and Stream.* files for convenience (Stream needed)
Ises=mis_read_ises(F.ext,F.path);
Blade=mis_read_blade(F.ext,F.path);
Stream=mis_read_stream(F.ext,F.path);

% Generates name and opens polarx file for reading
F.filename=['polarx.',F.ext];
F.file=fullfile(F.path,F.filename);
F.fid=fopen(F.file,'r');
if F.fid < 0; % Return error code if there is a problem
   warning(['Problem opening "',F.filename,'". Returning error code.']);
   Polarx=-1;
else % Read in the file
   % Read job title and and MISES version number
   t.head=fread(F.fid,2,'int16'); %fortran u/f header
   Polarx.title=fread(F.fid,32,'char')'; Polarx.title=char(Polarx.title); Polarx.title=deblank(Polarx.title);
   t.ddd=fread(F.fid,8,'char');
   Polarx.ver=fread(F.fid,1,'float64');
   t.foot=fread(F.fid,2,'int16'); %fortran u/f footer
   % Get number of Blades and number of Stream-wise points
   [Polarx.nbl,Polarx.npt]=mis_read_fortran(F.fid,2,'int32');
   % Get leading and trailing edge indicies
   t.ib=mis_read_fortran(F.fid,4.*Polarx.nbl,'int32');
   Polarx.ileb=t.ib(1:2:end-1);
   Polarx.iteb=t.ib(2:2:end);
   % Get number of aerofoil spline points and number of wake array points
   t.ibw=mis_read_fortran(F.fid,2.*Polarx.nbl,'int32');
   Polarx.nptb=t.ibw(1:2:end-1);
   Polarx.iwak=t.ibw(2:2:end);
   % Get blade spline coordinates
   Polarx.xb=cell(Polarx.nbl,1);
   Polarx.yb=cell(Polarx.nbl,1);
   for ii=1:Polarx.nbl
      t.xyb=mis_read_fortran(F.fid,2.*Polarx.nptb(ii),'float64');
      Polarx.xb{ii}=t.xyb(1:2:end-1);
      Polarx.yb{ii}=t.xyb(2:2:end);
   end
   % Read pressure distribution and BL params for each solution
   t.isoln=0; t.eof=0;
   while t.eof == 0 & t.isoln < F.nsmax
      t.isoln=t.isoln+1;
      t.head=fread(F.fid,2,'int16'); % fortran u/f header
      if ~isempty(t.head);
         % Read case summary
         t.nread=13 + 2.*Polarx.nbl;
		 raw.sinl(t.isoln)=fread(F.fid,1,'float64');
		 raw.sout(t.isoln)=fread(F.fid,1,'float64');
		 raw.minl(t.isoln)=fread(F.fid,1,'float64');
		 raw.mout(t.isoln)=fread(F.fid,1,'float64');
		 % dl467
% 		 raw.roinl(t.isoln)=fread(F.fid,1,'float64');
% 		 raw.qinl(t.isoln)=fread(F.fid,1,'float64');
% 		 raw.hinl(t.isoln)=fread(F.fid,1,'float64');
% 		 raw.pstinl(t.isoln)=fread(F.fid,1,'float64');
		 %
		 raw.p1pt(t.isoln)=fread(F.fid,1,'float64');
		 raw.p2pt(t.isoln)=fread(F.fid,1,'float64');
         raw.int(t.isoln)=fread(F.fid,1,'int32');
         raw.fturb(t.isoln)=fread(F.fid,1,'float64');
         raw.drvt(t.isoln)=fread(F.fid,1,'float64');
         raw.bvrn(1,t.isoln)=fread(F.fid,1,'float64');
         raw.bvrn(2,t.isoln)=fread(F.fid,1,'float64');
         raw.omega(t.isoln)=fread(F.fid,1,'float64');
         raw.omegv(t.isoln)=fread(F.fid,1,'float64');
         raw.xtr(:,t.isoln)=fread(F.fid,2.*Polarx.nbl,'float64');
         t.foot=fread(F.fid,2,'int16'); % fortran u/f footer
         % Read pressure distributions and boundary layer parameters
         t.nread=15.*Polarx.npt;
         for ii=1:2.*Polarx.nbl
            t.pbl=mis_read_fortran(F.fid,t.nread,'float64');
            t.pbl=reshape(t.pbl,15,Polarx.npt)';
            raw.x(:,ii,t.isoln)=t.pbl(:,1);
            raw.xi(:,ii,t.isoln)=t.pbl(:,2);
            raw.cp(:,ii,t.isoln)=t.pbl(:,3);
            raw.th(:,ii,t.isoln)=t.pbl(:,4);
            raw.dstr(:,ii,t.isoln)=t.pbl(:,5);
            raw.hbar(:,ii,t.isoln)=t.pbl(:,6);
            raw.uedg(:,ii,t.isoln)=t.pbl(:,7);
			raw.cf(:,ii,t.isoln)=t.pbl(:,8);
			raw.ctau(:,ii,t.isoln)=t.pbl(:,9);
			% new from Ho-on
			raw.hk(:,ii,t.isoln)=t.pbl(:,10);
			raw.hs(:,ii,t.isoln)=t.pbl(:,11);
			raw.rt(:,ii,t.isoln)=t.pbl(:,12);
			raw.di(:,ii,t.isoln)=t.pbl(:,13);
			raw.us(:,ii,t.isoln)=t.pbl(:,14);
			raw.cteq(:,ii,t.isoln)=t.pbl(:,15);
         end
      else
         t.eof=1;
         t.isoln=t.isoln-1;
      end
   end
   Polarx.nsoln=t.isoln;
   if Polarx.nsoln == 0 % Check for convergence failure
       Polarx=1;
       disp('Polar Convergence Problems');
       return
   end
   if t.isoln >= F.nsmax;
       warning('\n Max number of Polar solutions reached.  Change NSMAX to allow more solutions to be read if necessary.\n');
   end
   clear t
   % Sort and write files in ascending incidence
   [t.ddd,t.IX]=sort(raw.sinl);
   for ii=1:length(raw.sinl)
      Polarx.sinl(ii)=raw.sinl(t.IX(ii));
      Polarx.sout(ii)=raw.sout(t.IX(ii));
      Polarx.minl(ii)=raw.minl(t.IX(ii));
      Polarx.mout(ii)=raw.mout(t.IX(ii));
	  % dl467
% 	  if isfield(raw,'roinl') == 1 && isfield(raw,'qinl') == 1
% 		  Polarx.roinl(ii)=raw.roinl(t.IX(ii));
% 		  Polarx.qinl(ii)=raw.qinl(t.IX(ii));
% 		  Polarx.hinl(ii)=raw.hinl(t.IX(ii));
% 		  Polarx.pstinl(ii)=raw.pstinl(t.IX(ii));
% 	  end
	  %
      Polarx.p1pt(ii)=raw.p1pt(t.IX(ii));
      Polarx.p2pt(ii)=raw.p2pt(t.IX(ii));
      Polarx.int(ii)=raw.int(t.IX(ii));
      Polarx.fturb(ii)=raw.fturb(t.IX(ii));
      Polarx.drvt(ii)=raw.drvt(t.IX(ii));
      Polarx.bvrn(1,ii)=raw.bvrn(1,t.IX(ii));
      Polarx.bvrn(2,ii)=raw.bvrn(2,t.IX(ii));
      Polarx.omega(ii)=raw.omega(t.IX(ii));
      Polarx.omegv(ii)=raw.omegv(t.IX(ii));
      Polarx.xtr(:,ii)=raw.xtr(:,t.IX(ii));
      % Convert SINL and SOUT into angles measured in degrees
      Polarx.binl(ii)=atan(Polarx.sinl(ii)).*180./pi;
      Polarx.bout(ii)=atan(Polarx.sout(ii)).*180./pi;
      % Continue - these used to be 3D matricies, changed to cells by TH
      Polarx.x{ii}=raw.x(:,:,t.IX(ii));
      Polarx.xi{ii}=raw.xi(:,:,t.IX(ii));
      Polarx.cp{ii}=raw.cp(:,:,t.IX(ii));
      Polarx.th{ii}=raw.th(:,:,t.IX(ii));
      Polarx.dstr{ii}=raw.dstr(:,:,t.IX(ii));
      Polarx.hbar{ii}=raw.hbar(:,:,t.IX(ii));
      Polarx.uedg{ii}=raw.uedg(:,:,t.IX(ii));
      Polarx.cf{ii}=raw.cf(:,:,t.IX(ii));
      Polarx.ctau{ii}=raw.ctau(:,:,t.IX(ii));
	  Polarx.hk{ii}=raw.hk(:,:,t.IX(ii));
      Polarx.hs{ii}=raw.hs(:,:,t.IX(ii));
      Polarx.rt{ii}=raw.rt(:,:,t.IX(ii));
      Polarx.di{ii}=raw.di(:,:,t.IX(ii));
      Polarx.us{ii}=raw.us(:,:,t.IX(ii));
      Polarx.cteq{ii}=raw.cteq(:,:,t.IX(ii));
	  
	   % Calulate other BL parameters
      Polarx.cd{ii} = Polarx.di{ii}.*Polarx.hs{ii}/2; % Dissipation coefficient %%% DOESN'T WORK %%%
      Polarx.en{ii} = Polarx.hs{ii}.*Polarx.th{ii}; % Energy thickness
	  
      % Calculate surface length (interpolating for radius), see mises manual
      if isstruct(Stream); 
         Polarx.s{ii} = Polarx.xi{ii}.*interp1(Stream.x,Stream.r,Polarx.xi{ii},'spline','extrap');
      else; 
         error(' Need Stream file!'); 
      end;
      % Calculate xp and sp (note - uses coords from blade number one, but this should be ok)
      Polarx.xp{ii}=100.*((Polarx.x{ii}-min(Polarx.xb{1}))./(max(Polarx.xb{1})-min(Polarx.xb{1})));
      Polarx.sp{ii}(:,1)=100.*((Polarx.s{ii}(:,1)-Polarx.s{ii}(Polarx.ileb(1),1))./(Polarx.s{ii}(Polarx.iteb(1),1)-Polarx.s{ii}(Polarx.ileb(1),1)));
      Polarx.sp{ii}(:,2)=100.*((Polarx.s{ii}(:,2)-Polarx.s{ii}(Polarx.ileb(2),2))./(Polarx.s{ii}(Polarx.iteb(2),2)-Polarx.s{ii}(Polarx.ileb(2),2)));
      % Find ratio of specific heats
      if ii==1
         var.gam=1.4;
         t.dgam=1e-6;
         jj=0;
         while t.dgam >= 1e-6 & jj < 20;
            jj=jj+1;
            var.gamold=var.gam;
            var.fga=-log(1+(var.gam-1).*Polarx.minl(ii).*Polarx.minl(ii)./2)./log(Polarx.p1pt(ii));
            var.gam=1./(1-var.fga);
            t.dgam=abs(var.gam-var.gamold);
         end
      end
      % Convert read in values of Cp into isentropic surface Mach number
      var.cpref=1./(0.5.*var.gam.*Polarx.minl(ii).^2);
      var.p_pinl=Polarx.cp{ii}/var.cpref;
      var.p_pt=var.p_pinl.*Polarx.p1pt(ii);
      var.msq=2.*(var.p_pt.^(-var.fga)-1)./(var.gam-1);
      t.lt0=(var.msq<0);
      var.msq(t.lt0)=0;
      Polarx.mn{ii}=sqrt(var.msq);
      % Convert read in values of Cp into actual values of Cp
      Polarx.cp{ii}=Polarx.cp{ii}-var.cpref;
   end
   % Order struct
   Polarx=orderfields(Polarx);
   % Close idat.* file
   t.isclosed=fclose(F.fid);
   if t.isclosed < 0; warning(['Problem closing "',F.filename,'"']); end
   clear t
end

return