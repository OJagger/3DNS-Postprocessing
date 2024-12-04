function Idat=mis_read_idat(varargin)
% Function to read Mises idat.xxx files.
% --
% Initialise variables
switch nargin; % Deal with inputs
   case {0}; F.ext='mises'; F.path=pwd; % Assume current directory and file extension 'mises'
   case {1}; F.ext=varargin{1}; F.path=pwd;
   case {2}; F.ext=varargin{1}; F.path=varargin{2};
end
F.path = strrep(F.path,'TURBOSTREAM','MISES');
clear varargin
% fprintf(' Reading idat file...\n');

% Generate name and open idat file for reading
F.filename=['idat.',F.ext]; F.file=fullfile(F.path,F.filename); F.fid=fopen(F.file,'r');
if F.fid < 0; % Return error code if there is a problem
   warning(['Problem opening "',F.filename,'". Returning error code.']);
   Idat=[];
else % Read in the file
   % Read job title and and MISES version number
   t.head=fread(F.fid,2,'int16'); % fortran u/f header
   Idat.title=fread(F.fid,32,'char')'; Idat.title=char(Idat.title); Idat.title=deblank(Idat.title);
   t.foot=fread(F.fid,2,'int16'); % fortran u/f footer
   % Read state variables
   [Idat.nstati,Idat.nstatr] = mis_read_fortran(F.fid,2,'int32');
   Idat.istate = mis_read_fortran(F.fid,Idat.nstati,'int32');
   Idat.rstate = mis_read_fortran(F.fid,Idat.nstatr,'float64');
   % Unpack useful constants
   Idat.xinl = Idat.rstate(49);
   Idat.xout = Idat.rstate(50);
   Idat.pitch = Idat.rstate(51);
   Idat.sinl = Idat.rstate(23);
   Idat.binl = atan(Idat.sinl)*180/pi;
   % Unpack useful indicies
   Idat.ii = Idat.istate(1);
   Idat.jj = Idat.istate(2);
   Idat.nbl = Idat.istate(3);
   ns = Idat.nbl.*2;
   iih = Idat.istate(5);
   iip = Idat.istate(6);
   nbitn = Idat.istate(8);
   ngmode = Idat.istate(13);
   % Read mesh indicies
   store = mis_read_fortran(F.fid,4.*ns,'int32');
   store = reshape(store,4,ns);
   Idat.jbld = store(1,:);
   Idat.ninl = store(2,:);
   Idat.nbld = store(3,:);
   Idat.nout = store(4,:);
   store = mis_read_fortran(F.fid,3.*Idat.nbl,'int32');
   store = reshape(store,3,Idat.nbl);
   Idat.iib = store(1,:);
   Idat.ible = store(2,:);
   Idat.nwak = store(3,:);
   % Read blade coordinates and inlet/ outlet grid spacing arrays
   for ib = 1:Idat.nbl
      store = mis_read_fortran(F.fid,5.*Idat.iib(ib),'float64');
      store = reshape(store,5,Idat.iib(ib));
      Idat.xb{ib}  = store(1,:);
      Idat.yb{ib}  = store(2,:);
      Idat.xpb{ib} = store(3,:);
      Idat.ypb{ib} = store(4,:);
      Idat.sb{ib}  = store(5,:);
      store = mis_read_fortran(F.fid,5.*Idat.ii,'float64');
      store = reshape(store,5,Idat.ii);
      Idat.sginl{ib} = store(1,:);
      Idat.sgout{ib} = store(2,:);
      Idat.xw{ib}    = store(3,:);
      Idat.yw{ib}    = store(4,:);
      Idat.wgap{ib}  = store(5,:);
   end
   % Read streamline/streamtube/blade adjustment modes
   store = mis_read_fortran(F.fid,4,'int32');
   Idat.nbvrx = store(1);
   Idat.nmodx = store(2);
   Idat.nparx = store(3);
   Idat.nmovx = store(4);
   if Idat.nbvrx > 0; Idat.bvrn = mis_read_fortran(F.fid,Idat.nbvrx,'float64'); end
   if Idat.nmodx > 0; Idat.modn = mis_read_fortran(F.fid,Idat.nmodx,'float64'); end
   if Idat.nparx > 0; Idat.parn = mis_read_fortran(F.fid,Idat.nparx,'float64'); end
   if Idat.nmovx > 0; Idat.movn = mis_read_fortran(F.fid,Idat.nmovx,'float64'); end
   % Read grid coordinates and mass fraction between streamlines
   for jj = 1:Idat.jj
      % Mass fraction between streamlines
      Idat.mfrac(jj) = mis_read_fortran(F.fid,1,'float64');
      % Grid coordinates
      store = mis_read_fortran(F.fid,3.*Idat.ii,'float64');
      store = reshape(store,3,Idat.ii)';
      Idat.x(:,jj) = store(:,1);
      Idat.y(:,jj) = store(:,2);
      Idat.ro(:,jj) = store(:,3);
      
      %% dl467
      if jj ~= 1 
      Idat.x_cell(:,jj-1) = 0.5 * ( Idat.x(1:end-1,jj)+Idat.x(1:end-1,jj-1) );
      Idat.y_cell(:,jj-1) = 0.5 * ( Idat.y(1:end-1,jj)+Idat.y(1:end-1,jj-1) );
      end
      %%
   end
   
   for ii = 1:Idat.ii
      if ii ~= 1 
      Idat.x_cell(ii-1,:) = 0.5 * ( Idat.x(ii,1:end-1)+Idat.x(ii-1,1:end-1) );
      Idat.y_cell(ii-1,:) = 0.5 * ( Idat.y(ii,1:end-1)+Idat.y(ii-1,1:end-1) );
      end
   end
   
    %%
   % Calculate xp for each surface (percent chord on blade surfaces), based on
   % the geometry of blade 1
   Idat.xp(:,1)=100.*((Idat.x(:,1))./(max(Idat.xb{1})-min(Idat.xb{1})));
   Idat.xp(:,2)=100.*((Idat.x(:,end))./(max(Idat.xb{1})-min(Idat.xb{1})));
   % Read boundary layer and wall parameters
   for is=1:ns
      store=mis_read_fortran(F.fid,5.*Idat.ii,'float64');
      store=reshape(store,5,Idat.ii)';
      Idat.sg(:,is)=store(:,1);
      Idat.disp(:,is)=store(:,2);
      Idat.sdisp(:,is)=store(:,3);
      Idat.pspec(:,is)=store(:,4);
      Idat.dhsdhb(:,is)=store(:,5);
      store=mis_read_fortran(F.fid,5.*Idat.ii,'float64');
      store=reshape(store,5,Idat.ii)';
      Idat.th(:,is)=store(:,1);
      Idat.dstr(:,is)=store(:,2);
      Idat.uedg(:,is)=store(:,3);
      Idat.ctau(:,is)=store(:,4);
      Idat.tau(:,is)=store(:,5);
      store=mis_read_fortran(F.fid,4.*Idat.ii,'float64');
      store=reshape(store,4,Idat.ii)';
      Idat.hwall(:,is)=store(:,1);
      Idat.qwall(:,is)=store(:,2);
      Idat.mwall(:,is)=store(:,3);
      Idat.psiw(:,is)=store(:,4);
   end
   % Read meridional and streamtube details
   npl=(3 + 2.*(Idat.nbvrx+1)); % number of elements "per line"
   store=mis_read_fortran(F.fid,npl.*iih,'float64');
   store=reshape(store,npl,iih)';
   Idat.xh=store(:,1);
   Idat.rh=store(:,2);
   Idat.rph=store(:,3);
   iist=4;
   iend=iist+Idat.nbvrx;
   Idat.bh=store(:,iist:iend);
   iist=iend+1;
   iend=iist+Idat.nbvrx;
   Idat.bph=store(:,iist:iend);
   Idat.ibhdef=mis_read_fortran(F.fid,Idat.nbvrx,'int32');
   % Read specified pressure loss
   store=mis_read_fortran(F.fid,3.*iip,'float64');
   store=reshape(store,3,iip)';
   Idat.xrl=store(:,1);
   Idat.prl=store(:,2);
   Idat.prlx=store(:,3);
   %Read forces and moments
   store=mis_read_fortran(F.fid,5.*Idat.nbl,'float64');
   store=reshape(store,5,Idat.nbl)';
   Idat.bldfx=store(:,1);
   Idat.bldfy=store(:,2);
   Idat.bldmz=store(:,3);
   Idat.bldvl=store(:,4);
   Idat.ppsgw=store(:,5);
   % Get spline properties
   store=mis_read_fortran(F.fid,2.*Idat.nbl,'float64');
   store=reshape(store,2,Idat.nbl)';
   Idat.sble=store(:,1);
   Idat.sblold=store(:,2);
   store=mis_read_fortran(F.fid,5.*Idat.nbl,'float64');
   store=reshape(store,5,Idat.nbl)';
   Idat.sblegn=store(:,1);
   Idat.xblegn=store(:,2);
   Idat.yblegn=store(:,3);
   Idat.sstg=store(:,4);
   Idat.swak=store(:,5);
   % Read transition and 2nd derivatives of pressures (left hand and right hand respectivley)
   store = mis_read_fortran(F.fid,3.*is,'float64');
   store = reshape(store,3,is)';
   Idat.xtr  = store(:,1);
   Idat.pxx0 = store(:,2);
   Idat.pxx1 = store(:,3);
   store = mis_read_fortran(F.fid,2.*is,'int32');
   store = reshape(store,2,is)';
   Idat.itran = store(:,1);
   Idat.ktran = store(:,2);
   % Read suction details
   Idat.nsuct = mis_read_fortran(F.fid,1,'int32');
   if Idat.nsuct > 0
      Idat.issuct = mis_read_fortran(F.fid,Idat.nsuct,'int32');
      store = mis_read_fortran(F.fid,3.*Idat.nsuct,'float32');
      store = reshape(store,3,Idat.isuct);
      Idat.cqsuct = store(1,:);
      Idat.sgsuct(1,:) = store(2,:);
      Idat.sgsuct(2,:) = store(3,:);
   end
   % Read packed isentropic-cell bits
   for ib = 1:nbitn
      Idat.isbits(:,ib) = mis_read_fortran(F.fid,Idat.ii,'int32');
   end
   % Read mode specification lines
   Idat.ngmode = mis_read_fortran(F.fid,1,'int32');
   if Idat.ngmode > 0
      store = mis_read_fortran(F.fid,3.*ngmode);
      Idat.igmode = reshape(store,3,ngmode);
      store = mis_read_fortran(F.fid,3.*ngmode);
      Idat.agmode = reshape(store,3,ngmode);
   end
   % Read convergence and choke flags
   [Idat.lconv,Idat.lchoke] = mis_read_fortran(F.fid,2,'int32');
   
   % Read loss and turning
   try
       [Idat.omeg,Idat.sout] = mis_read_fortran(F.fid,2,'float64');
   catch
   end
   
   %% dl467 determine Mach number, loss and other primary variable contours
   ga = 1.4;
   Idat.Mis = sqrt(2/(ga-1)*(Idat.ro(1:end-1,1:end-1).^-(ga-1)-1));
   % determine area normal to streamline for each cell
   Idat.dx_i_line = double(diff(Idat.x,1,2));
   Idat.dy_i_line = double(diff(Idat.y,1,2));
   Idat.dx_j_line = double(diff(Idat.x,1,1));
   Idat.dy_j_line = double(diff(Idat.y,1,1));
   Idat.di_line(:,:,1) =  Idat.dx_i_line(1:end-1,:);
   Idat.di_line(:,:,2) =  Idat.dy_i_line(1:end-1,:);
   Idat.dj_line(:,:,1) =  Idat.dx_j_line(:,1:end-1);
   Idat.dj_line(:,:,2) =  Idat.dy_j_line(:,1:end-1);
   
    
   % calculate unit normal vector to streamwise distance
   Idat.es(:,:,1) = double((Idat.dj_line(:,:,1)./sqrt(sum(Idat.dj_line.^2,3))));
   Idat.es(:,:,2) = double((Idat.dj_line(:,:,2)./sqrt(sum(Idat.dj_line.^2,3))));
   Idat.en(:,:,1) = -Idat.es(:,:,2);
   Idat.en(:,:,2) = Idat.es(:,:,1);
   Idat.dn_line(:,:,1) = double(double(dot(Idat.di_line,Idat.en,3)).*Idat.en(:,:,1));
   Idat.dn_line(:,:,2) = double(double(dot(Idat.di_line,Idat.en,3)).*Idat.en(:,:,2));
   Idat.dn_line_norm(:,:,1) = Idat.dn_line(:,:,1)./sqrt(sum(Idat.dn_line.^2,3));
   Idat.dn_line_norm(:,:,2) = Idat.dn_line(:,:,2)./sqrt(sum(Idat.dn_line.^2,3));
   % calculate normal area to streamline of each cell
   % take radial contraction into account
   [Idat.xh,i_uniq] = unique(Idat.xh);
   Idat.Ax_ratio = interp1(Idat.xh,Idat.bh(i_uniq),Idat.x_cell);
   Idat.dAn = double(sqrt(sum(Idat.dn_line.^2,3)).*Idat.Ax_ratio);
   
   % determine inlet area and hence velocity fraction
   Idat.An = double(sum(Idat.dAn,2));
   Idat.An1 = double(Idat.An(1));
   Idat.ro1 = double(sum(Idat.ro(1,1:end-1))/size(Idat.ro(:,1:end-1),2)); 
   for j = 1:size(Idat.dAn,2)
   Idat.Vfrac(:,j) = double(Idat.mfrac(j)./(Idat.dAn(:,j)./Idat.An1).*Idat.ro1./Idat.ro(1:end-1,j));
   end
   %% CAUTION: assumption no change in radius for Torel to be constant (V_cpTo = f(M))
   Idat.Mfrac = double(Idat.Vfrac./sqrt(1+(ga-1)/2*(1-Idat.Vfrac.^2)));
   %%
   Idat.Minl = double(sum(Idat.Mis(1,1:end-1))./(size(Idat.Mis,2)-1));
   Idat.M = double(Idat.Mfrac.*Idat.Minl);
   Idat.Po_ratio = double(( (1+(ga-1)/2.*Idat.M.^2)./(1+(ga-1)/2.*Idat.Mis.^2) ).^(ga/(ga-1)));
   Idat.dPo_Po = double(Idat.Po_ratio - 1);
   %%
   
   % Order struct
   Idat=orderfields(Idat);
   % Close idat.* file
   t.isclosed=fclose(F.fid);
   if t.isclosed < 0; warning(['Problem closing "',F.filename,'"']); end
   clear t
end

return