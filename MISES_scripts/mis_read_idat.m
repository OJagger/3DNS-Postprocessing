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
      Idat.ro_cell(:,jj-1) = 0.5* ( Idat.ro(1:end-1,jj)+Idat.ro(1:end-1,jj-1) );
      end
      %%
   end
   
   for ii = 1:Idat.ii
      if ii ~= 1 
      Idat.x_cell(ii-1,:) = 0.5 * ( Idat.x(ii,1:end-1)+Idat.x(ii-1,1:end-1) );
      Idat.y_cell(ii-1,:) = 0.5 * ( Idat.y(ii,1:end-1)+Idat.y(ii-1,1:end-1) );
      Idat.ro_cell(ii-1,:) = 0.5* ( Idat.ro(ii,1:end-1)+Idat.ro(ii-1,1:end-1) );
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
       disp('')
   end
   
   %% dl467 determine Mach number, loss and other primary variable contours
   ga = 1.394;
%    Idat.Mis = sqrt(2/(ga-1)*(Idat.ro_cell.^-(ga-1)-1));
   Idat.Mis = sqrt(2/(ga-1)*(Idat.ro(1:end-1,1:end-1).^(-(ga-1))-1));
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
   Idat.dAn_Ax1_dl = double(sqrt(sum(Idat.dn_line.^2,3)));
   Idat.dAn = double(sqrt(sum(Idat.dn_line.^2,3)).*Idat.Ax_ratio);
   Idat.dAn_A1 = Idat.dAn./repmat(Idat.dAn(1,:),size(Idat.dAn,1),1);
   
   % based on Drela thesis definitions found in pg. 22-23
   Ax1 = double(0.5*(Idat.dx_i_line(1:end-1,:)+Idat.dx_i_line(2:end,:)));
   Ay1 = double(0.5*(Idat.dy_i_line(1:end-1,:)+Idat.dy_i_line(2:end,:)));
   Ax2 = double(0.5*(Idat.dx_i_line(2:end-1,:)+Idat.dx_i_line(3:end,:)));
   Ay2 = double(0.5*(Idat.dy_i_line(2:end-1,:)+Idat.dy_i_line(3:end,:)));
   sx1 = double(0.5*(Idat.dx_j_line(:,2:end)+Idat.dx_j_line(:,1:end-1)));
   sy1 = double(0.5*(Idat.dy_j_line(:,2:end)+Idat.dy_j_line(:,1:end-1)));
   sx2 = double(0.5*(Idat.dx_j_line(2:end,2:end)+Idat.dx_j_line(2:end,1:end-1)));
   sy2 = double(0.5*(Idat.dy_j_line(2:end,2:end)+Idat.dy_j_line(2:end,1:end-1)));
   s1 = double(sqrt(sx1.^2+sy1.^2));
   s2 = double(sqrt(sx2.^2+sy2.^2));
   Idat.dAn_Ax1 = abs(double((sx1.*Ay1-sy1.*Ax1)./s1));
   Idat.dAn_Ax2 = abs(double((sx2.*Ay2-sy2.*Ax2)./s2));
   Idat.dAn = Idat.dAn_Ax1.*Idat.Ax_ratio;
   
   % determine inlet area and important inlet quantities
   Idat.An = double(sum(Idat.dAn,2));
   Idat.An1 = double(Idat.An(1));
   Idat.ro1_area = double(sum(Idat.ro(1,1:end-1))/size(Idat.ro(:,1:end-1),2)); % stagnation density is equal to unity so this is the true density value
   Idat.ro_stag = 1; % MISES definition
   %% NOTE THIS IS FOR ZERO CHANGE IN RADIUS
   Idat.h_stag = 1/(ga-1); % MISES definition for zero change in radius
   Idat.cp = 1003.308;
   Idat.To_stag = Idat.h_stag./Idat.cp;
   %%
   Idat.q1_area = double(sqrt(2.*Idat.h_stag.*(1-(Idat.ro1_area./Idat.ro_stag).^(ga-1)))); % eq.2.29 pg. 31
   Idat.q1 = double(sqrt(2.*Idat.h_stag.*(1-(Idat.ro(1,1:end-1)./Idat.ro_stag).^(ga-1)))); % eq.2.29 pg. 31
   Idat.Minl = double(sum(Idat.Mis(1,1:end-1))./(size(Idat.Mis,2)-1));
%    Idat.m_tot = Idat.An1.*Idat.ro1.*Idat.V1;
   Idat.dm_tot = Idat.dAn(1,:).*Idat.ro(1,1:end-1).*Idat.q1;
   Idat.m_tot = sum(Idat.dm_tot);
   
   % calcualte strematube parameters which can be plotted out from MISES
   Idat.dAn_A1 = Idat.dAn./repmat(Idat.dAn(1,:),size(Idat.dAn,1),1); % normal area ratio
   Idat.ro_ro1 = Idat.ro(1:end-1,1:end-1)./repmat(Idat.ro(1,1:end-1),size(Idat.ro,1)-1,1); % density ratio
   Idat.q_q1 = 1./(Idat.dAn_A1.*Idat.ro_ro1); % velocity ratio
   Idat.q = Idat.m_tot.*repmat(Idat.mfrac(1:end-1),size(Idat.ro,1)-1,1)./(Idat.dAn.*Idat.ro(1:end-1,1:end-1));
   Idat.P = (ga-1)./ga.*Idat.ro(1:end-1,1:end-1).*(Idat.h_stag-0.5*(Idat.q_q1.*repmat(Idat.q1,size(Idat.q_q1,1),1)).^2); % eq. 2.12 pg 25
   Idat.P1 = double(sum(Idat.P(1,:))/size(Idat.P,2));
%    Idat.Po_stag = Idat.P1.*(1+(ga-1)/2.*Idat.Minl.^2).^(ga/(ga-1)); %
%    this is a check
   Idat.Po_stag = 1./ga; % from mises definition
   Idat.P_Po1 = Idat.P./Idat.Po_stag;
   
   % determine velocity components
   Idat.Vx = Idat.q.*sx1./s1;
   Idat.Vt = Idat.q.*sy1./s1;

   % claculate true Mach number
   Idat.M = double(sqrt(Idat.ro(1:end-1,1:end-1).*Idat.q.^2./(ga.*Idat.P))); % M^2=rhou^2/gaP
   % determine total pressure loss
   Idat.Po_ratio = double(( (1+(ga-1)/2.*Idat.M.^2)./(1+(ga-1)/2.*Idat.Mis.^2) ).^(ga/(ga-1)));
   Idat.dPo_Po = double(Idat.Po_ratio - 1);
   % correct beacuse freestream loss is sometimes predicted above the
   % shock - this is because of precision error at the two ends
   Idat.dPo_Po_corr = Idat.dPo_Po;
   Idat.dPo_Po_corr(Idat.dPo_Po_corr>0) = 0;
   Idat.Po = Idat.Po_ratio.*Idat.Po_stag;
   Idat.To = ones(size(Idat.Po_ratio,1),size(Idat.Po_ratio,2)).*Idat.To_stag;
   
   
   % calculate supersonic solutions to area change
   % Mass flow function: m * (cp * To)^0.5 / (A * Po)
   % Range of Mach numbers
   Mach = 0:0.0005:2;
   M_func = (ga / ((ga-1)^0.5)) * Mach .* ((1 + 0.5 * (ga-1) * Mach.^2) .^ (-0.5 * (ga+1) / (ga-1)));
   
   % Extract separate subsonic and supersonic mass flow function solutions
   q_sub = Mach <= 1; Mach_sub = Mach(q_sub); M_func_sub = M_func(q_sub);
   q_sup = Mach >= 1; Mach_sup= Mach(q_sup); M_func_sup = M_func(q_sup);
   
   % calculate isentropic area change from eq. from Bernstein et. al.
   Idat.dAn_is = double( 1./(ga./sqrt(ga-1).*repmat(Idat.Po_stag./(Idat.dm_tot.*sqrt(Idat.h_stag)),size(Idat.ro,1)-1,1).*Idat.Mis.*(1+(ga-1)/2.*Idat.Mis.^2).^(-0.5.*(ga+1)./(ga-1))) );
%    Idat.dAn_is = Idat.dAn_A1_is.*repmat(Idat.dAn(1,:),size(Idat.dAn,1),1);
   Idat.An_is = double(sum(Idat.dAn_is,2));
   
   
   % calcualte superosnic isentropic exit solution for the same change in
   % individual streamtube area
   m_nondim  = double(repmat(Idat.dm_tot.*sqrt(Idat.h_stag),size(Idat.dAn,1),1)./(Idat.dAn_is.*Idat.Po_stag));
   for j=1:size(m_nondim,2)
       [~,i_min]=min(abs(repmat(m_nondim(:,j),[1 size(M_func_sup,2)])'-repmat(M_func_sup,[size(m_nondim,1) 1])'));
       Idat.Mis_ss(:,j)=Mach_sup(i_min);
       [~,i_max(j)]=max(Idat.Mis(:,j));
       Idat.Mis_ss_temp = Idat.Mis_ss(i_max(j)+1:end,j);
       Idat.Mis_temp = Idat.Mis(i_max(j)+1:end,j);
       Idat.Mis_ss_exit(:,j) = [Idat.Mis(1:i_max(j),j) ;  Idat.Mis_temp(Idat.Mis_temp>=1); Idat.Mis_ss_temp(Idat.Mis_temp<1)];
   end

   % streamwise interpolate necessary parameters to calculate the beta
   % parameter
   if isnan(Idat.sinl) == 0 % check if solution has not run
   Idat.s = [zeros(1,size(Idat.x,2)) ; cumsum((diff(Idat.x,1,1).^2+diff(Idat.y,1,1).^2).^0.5,1)];
   Idat.s = Idat.s(1:end-1,1:end-1); % stremawise distance
   Idat.s_datum = Idat.s(:,round(size(Idat.s,2)/2)); % this is used as the reference streamwise distance
   for j=1:size(Idat.s,2)
       Idat.dAn_interp(:,j) = interp1(Idat.s(:,j),Idat.dAn(:,j),Idat.s_datum);
       Idat.dAn_is_interp(:,j) = interp1(Idat.s(:,j),Idat.dAn(:,j),Idat.s_datum);
       Idat.Mis_interp(:,j) = interp1(Idat.s(:,j),Idat.Mis(:,j),Idat.s_datum);
       Idat.M_interp(:,j) = interp1(Idat.s(:,j),Idat.M(:,j),Idat.s_datum);
       Idat.Mis_ss_exit_interp(:,j) = interp1(Idat.s(:,j),Idat.Mis_ss_exit(:,j),Idat.s_datum);
       Idat.P_interp(:,j) = interp1(Idat.s(:,j),Idat.P(:,j),Idat.s_datum);
   end
   Idat.An_interp = double(sum(Idat.dAn_interp,2));
   Idat.An_is_interp = double(sum(Idat.dAn_is_interp,2));
    
   % calculate the beta factor from Compressible compound nozzles Bernstein
   % et.al.
   Idat.beta = double(sum(Idat.dAn_interp./ga.*(1./Idat.M_interp.^2-1),2));
   Idat.beta_is = double(sum(Idat.dAn_interp./ga.*(1./Idat.Mis_interp.^2-1),2));
   Idat.beta_ss = double(sum(Idat.dAn_is_interp./ga.*(1./Idat.Mis_ss_exit_interp.^2-1),2));
   % calculate what it would be if it were nozzle
   Idat.P_mass = double(sum(repmat(Idat.dm_tot,size(Idat.P_interp,1),1).*Idat.P_interp,2)./Idat.m_tot);
   Idat.Mis_mass = double(sum(repmat(Idat.dm_tot,size(Idat.Mis_interp,1),1).*Idat.Mis_interp,2)./Idat.m_tot);
   Idat.beta_1D = double(Idat.An_interp./ga.*(1./Idat.Mis_mass.^2-1));
   
   % calculate pressure gradient in s and n directions
   % for j=1:size(Idat.s,2)
   %     [Idat.dP_ds(:,j),~] =  grad_mg(Idat.s_datum,Idat.P_interp(:,j));
   % end
   % for i=1:size(Idat.s,1)
   %     [Idat.dP_dn(i,:),~] =  grad_mg(cumsum(Idat.dAn_interp(i,:)),Idat.P_interp(i,:));
   % end

   % calculate outlet mass averaged quantities
   Idat.P_mass = double(sum(repmat(Idat.dm_tot,size(Idat.P_interp,1),1).*Idat.P_interp,2)./Idat.m_tot);
   Idat.P_area = double(sum(Idat.dAn_interp.*Idat.P_interp,2)./Idat.An_interp);
   Idat.dPo2_Po_mass = double(sum(Idat.dm_tot.*Idat.dPo_Po_corr(end,:))./Idat.m_tot);
   Idat.dPo2_Po_mass_datum = double(sum(Idat.dm_tot.*Idat.dPo_Po(end,:))./Idat.m_tot);
   Idat.Po2_mass_datum = double(sum(Idat.dm_tot.*Idat.Po(end,:))./Idat.m_tot);
   end
   %%
   
   % Order struct
   Idat=orderfields(Idat);
   % Close idat.* file
   t.isclosed=fclose(F.fid);
   if t.isclosed < 0; warning(['Problem closing "',F.filename,'"']); end
   clear t
end

return