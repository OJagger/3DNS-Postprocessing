function [p,h] = mis_plot_section(directory,h,col,plot_stuff,job,inc)
% MIS_PLOT_SECTION  Plot a range of data from a MISES run on a blade section
%
%   [p,h] = MIS_PLOT_SECTION(directory,h,col,plot_stuff,t)
%
%   directory - string of output file directory
%   h - optional figure handle
%   col - optional colour RGB vectors
%   plot_stuff - 0 or 1 fgeditor showing working
%   t - optional blade section information
%   p - output data structure

% Default to show plot
if exist('plot_stuff','var') == 0
    plot_stuff = 1;
end

% Default colour to black
if exist('col','var') == 0 || isempty(col) == 1
    col = [0 0 0];
end

directory = strrep(directory,'TURBOSTREAM','MISES');

%% Measure flow features from converged MISES solution

% Read in section parameters if not specified
% if exist('c','var') == 0 
load([directory 'section.mat']);
% end

% Read in flow file
if exist([directory 'polarx.mises'],'file') ~= 0
    [Polarx, Ises] = mis_read_polarx('mises',directory);
else
    disp('File Not Found')
    p = [];
    return
end

% Check the point is converged
if isfield(Polarx,'binl') == 0
    disp('Run Not Converged')
    p = [];
    return
end


% Find the design incidence from the polar run
if exist('inc','var')==0 || isempty(inc) == 1
	[~,j] = min(abs(Polarx.binl-Ises.binl));
else
    [~,j] = min(abs(Polarx.binl-(Ises.binl+inc)));
    disp(['Actual incidence is ' num2str(Polarx.binl(j)-Ises.binl) ])
end

% Record the outlet flow angle
p.Alpha = atand(Polarx.sout(j));

% Find the closest point to peak suction
Cp = Polarx.cp{j}(:,1); s = Polarx.s{j}(:,1) / Polarx.s{j}(Polarx.iteb(1),1);
Mis = Polarx.mn{j}(:,1);
q = s > 0.15 & s < 0.99; Cp = Cp(q); s = s(q); Mis = Mis(q);
[~, i] = min(Cp); i = max(i,3);


% Fit a polynomial through points near peak suction to find the maximum more accurately
s_temp = linspace(s(i-2),s(i+2),1000);
Q = polyval(polyfit(s(i-2:i+2),Cp(i-2:i+2),3),s_temp);
[~,i_max] = min(Q); p.s_Cp_max= s_temp(i_max);
M = polyval(polyfit(s(i-2:i+2),Mis(i-2:i+2),3),s_temp);
[~,i_max] = min(M); p.s_M_max= s_temp(i_max);
p.M_max = M(i_max);


% Get all idat filenames
A = [dir([directory 'idat.mises*'])];
F = cell(length(A),1); for n = 1:length(A); F{n} = A(n).name; end;

if exist('inc','var')==0 || isempty(inc) == 1
    % Read in the grid coodinates
    if exist([directory 'idat.mises_01'],'file') ~= 0
        Idat = mis_read_idat('mises_01',directory);
    else
        Idat = mis_read_idat('mises',directory);
    end
else
    for i=1:length(F)
        Idat = mis_read_idat(F{i}(6:end),directory);
        alpha_inlet(i) = atand(Idat.sinl);
    end
    [~,i_inc] = min(abs(alpha_inlet-(Ises.binl+inc)));
    Idat = mis_read_idat(F{i_inc}(6:end),directory);
end

% Check the point is converged
if isnan(Idat.binl) == 1
    disp('Run Not Converged')
    p = [];
    return
end

% extract bl stream
bl_stream = [Idat.x(Idat.ninl(1):end,1) Idat.y(Idat.ninl(1):end,1)].*job.r_mid;


% Calculate the relative angle from the leading edge centre
mt_stag = [Idat.x(Idat.ninl(2),end) Idat.y(Idat.ninl(2),end)];
phi_stag = atan2(c.mt_le_cen(2) - mt_stag(2),c.mt_le_cen(1) - mt_stag(1)) * 360 / (2*pi);
local_inc = phi_stag - c.chi_le;
p.local_inc = local_inc;
mt = c.mt(c.mt(:,1) < c.mt_le_cen(1),:);
[~,i_ss]=min(c.mt(:,1)-c.mt_le(1));
% c.mt_ss = flip(c.mt(1:i_ss,:));
c.mt_ss = c.mt(1:i_ss,:);
c.mt_ps = c.mt(i_ss+1:end,:);

mt_ss = c.mt(c.mt_ss(:,1) > c.mt_le_cen(1),:);
mt_ss = flip(mt_ss);
mt_ps = c.mt(c.mt_ps(:,1) > c.mt_le_cen(1),:);
mt_ps = flip(mt_ps);
% figure; hold on;
% plot(mt_ss(:,1),mt_ss(:,2),'-r')
% plot(mt_ps(:,1),mt_ps(:,2),'-b')
% plot(mt_ps(1,1),mt_ps(1,2),'xg')
% plot(mt_ss(:,1),mt_ss(:,2),'-r')

scl_ss = [0 ; cumsum(sum(diff(mt_ss,1,1).^2,2).^0.5,1)];
scl_ss = scl_ss/max(scl_ss);
scl_ps = [0 ; cumsum(sum(diff(mt_ps,1,1).^2,2).^0.5,1)];
scl_ps = scl_ps/max(scl_ps);
phi = atan2(c.mt_le_cen(2) - mt(:,2),c.mt_le_cen(1) - mt(:,1)) * 360 / (2*pi);


% Calculate the angle the stagnation point makes with the leading edge
[phi,i] = unique(phi); mt = mt(i,:);
dtdx = grad_mg(mt(:,1),mt(:,2));
psi = atand(dtdx);
dtdx_ss = grad_mg(mt_ss(:,1),mt_ss(:,2));
psi_ss = atand(dtdx_ss);
dtdx_ps = grad_mg(mt_ps(:,1),mt_ps(:,2));
psi_ps = atand(dtdx_ps);


% Correct the angle by the leading edge metal angle
p.psi_stag = 90 - (c.chi_le - interp1(phi,psi,phi_stag,'pchip'));
psi_stag = interp1(phi,psi,phi_stag,'pchip');

% Record the boundary layer parameters
q_1 = Polarx.ileb(1):Polarx.iteb(1);
q_2 = Polarx.ileb(2):Polarx.iteb(2);

s_1 = Polarx.s{j}(q_1,1); 
s_1 = s_1 / max(s_1);
s_2 = Polarx.s{j}(q_2,2); s_2 = s_2 / max(s_2);

Cp_1 = Polarx.cp{j}(q_1,1); Cp_2 = Polarx.cp{j}(q_2,2);
Mis_1 = Polarx.mn{j}(q_1,1); Mis_2 = Polarx.mn{j}(q_2,2);
Hb_1 = Polarx.hbar{j}(q_1,1);
Th_1 = Polarx.th{j}(q_1,1); Th_2 = Polarx.th{j}(q_2,2); 

% Record the shape factor profile parameters
p.s = s_1; p.Hb = Hb_1; p.q_1 = q_1;
p.Mis_1 = Mis_1;
p.M_max_ss = max(Mis_1);

%% calculate AtA1
if exist('job','var') == 1
% 	g=ts_read_hdf5([job.rjm_directory job.outname]);
% 	nj_mid = round(0.5*g{1}.attribute.nj);
	pitch = 2*pi*job.r_mid/c.N;
    %% legact chagne
    if isfield(job,'tchord')==0
        job.tchord = job.chord;
    end
    %%
	c.s_c = round(100*pitch/job.tchord)/100;
    %% legact chagne
    if isfield(job,'thick_max')==0
        c.t_c = job.t.thick_chord;
    else
        c.t_c = job.thick_max./job.tchord;
    end
    
    if plot_stuff == 1 || plot_stuff == 2
    
    if isfield(job,'dev') == 1
        % construct section WITH trailing edge
        ote = 0; dev =[]; plot_stuff = 0;
        cam_spline = 1;
        [xrt,xrt_cam,xrt_ps,xrt_ss]=bl_construct_section(job,plot_stuff,ote,dev,cam_spline);
        %% removed because it takes a long time
        %% calculate value of o_scosa
        [~,o_s]=ts_throat_calc(xrt_ss,xrt_ps,pitch,[],0);
        if exist('inc','var')==0 || isempty(inc) == 1
            c.oscosa=o_s/cosd(job.alpha_rel_inlet);
        else
            c.oscosa=o_s/cosd(abs(job.alpha_rel_inlet)+inc);
        end
        %%
        [rad_contr] = mis_radial_contraction(directory);
        c.AtA1 = c.oscosa*rad_contr;
    else
        % construct section WITH trailing edge
        ote = 0; dev =[];
        cam_spline = 1;
        %% this is legacy
        if isfield(job,'rad_le')==0
            [xrt,sl_cx,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
            xrt_ss = mca.xrt_ss; xrt_ps = mca.xrt_ps;
            %%
        else
            [xrt,xrt_cam,xrt_ps,xrt_ss]=bl_construct_section(job,0,ote,dev,cam_spline);
        end
        % calculate value of o_scosa
        bl = []; pitch_ref = pitch; plot_throat = 0; mode = directory; % use normal to midpassage streamline to get throat 
        alpha_inlet = []; xrt_stag_line = [];
%         [~,o_s,xrt_throat]=ts_throat_calc(xrt_ss,xrt_ps,pitch,bl,plot_throat,pitch_ref,xrt_stag_line,alpha_inlet,mode);
%         bl_stream(:,2) =  bl_stream(:,2) + pitch;
%         [~,o_s_bl,xrt_throat_bl]=ts_throat_calc(bl_stream,xrt_ps,pitch,bl,plot_throat,pitch_ref,xrt_stag_line,alpha_inlet,mode);
    [rad_contr,AtA1,AtA1_roV,rad_contr_o,o_s,xrt_throat] = mis_radial_contraction(directory);
    bl_inc = 1;
    [rad_contr_bl,AtA1_bl,AtA1_roV_bl,rad_contr_o_bl,o_s_bl,xrt_throat_bl] = mis_radial_contraction(directory,0,bl_inc);
    
    
    if exist('inc','var')==0 || isempty(inc) == 1
        c.oscosa=o_s/cosd(job.alpha_rel_inlet);
        c.oscosa_bl=o_s_bl/cosd(job.alpha_rel_inlet);
    else
        c.oscosa=o_s/cosd(abs(job.alpha_rel_inlet)+inc);
        c.oscosa_bl=o_s_bl/cosd(job.alpha_rel_inlet+inc);
    end
    
    
    
	c.AtA1 = c.oscosa*rad_contr;
    c.AtA1_v2 = AtA1;
    c.AtA1_rhoV  = AtA1_roV;
    c.AtA1_bl = c.oscosa_bl*rad_contr_bl;
    c.AtA1_v2_bl = AtA1_bl;
    c.AtA1_rhoV_bl  = AtA1_roV_bl;
% 	c.AtA1 = c.oscosa;

	p.AtA1 = c.AtA1;
    p.AtA1_v2 = c.AtA1_v2;
    p.AtA1_rhoV  = c.AtA1_rhoV;
    p.AtA1_bl = c.AtA1_bl;
    p.AtA1_v2_bl = c.AtA1_v2;
    p.AtA1_rhoV_bl  = c.AtA1_rhoV_bl;
	
    end
    
    end
    
    
end


% dl467 (code taken from H_calc area)
% determine where linear H-bar should start
% in case minimum value lies somewhere around the LE
%% this is legacy
if isfield(c,'ss_sb')==0
   c.ss_sb = mca.ss_sb ;
end
%%
ss_sb = c.ss_sb;
mt_max = p.s_Cp_max;
[~,i_min]=min(Hb_1);
s_min=s_1(i_min);
if s_min<0.1 || s_min>0.8
    s_min=ss_sb;
    H_min = interp1(s_1,Hb_1,s_min,'spline');
end
% in the case we have a shock in the system to control H-bar after the
% shock
if max(Mis_1)>1
%     [H_min,i_min2]=min(H_fit(s>x_max+0.25)); % find location of second dip
    [H_min,i_min2]=min(Hb_1(s_1>mt_max+0.25)); % find location of second dip
    H_interp_temp = Hb_1(s_1>mt_max+0.25);
    s_temp = s_1(s_1>mt_max+0.25);
    if i_min2 ~= length(s_1(s_1>mt_max+0.25));
    s_min = s_temp(i_min2);
    else
        i_min2 = i_min2 - 5;
		if i_min2<1
			i_min2 = 1;
		end
        s_min = s_temp(i_min2);
        H_min = H_interp_temp(i_min2);
    end
end

if isempty(s_min) == 1
	s_min = 0.97;
elseif s_min>0.975
	s_min = 0.97;
end

% if isempty(H_min) == 1
% 	H_min = H_98;
% end


s_bl_min=s_min;

% Interpolate shape factor at increased resolution
%s_bl_min = 0.45; 
s_bl_max = 1 - 0.5 * c.thick_te * c.thick_max / c.tchord - 0.01; 
% s_bl_mid = 0.7;
s_bl_mid = 0.5*(s_bl_min+s_bl_max);
s_bl = linspace(s_bl_min,s_bl_max,100)'; Hb_bl = interp1(s_1,Hb_1,s_bl,'pchip');
Th_bl = interp1(s_1,Th_1,s_bl,'pchip');

p.Th_bl = Th_bl;

% Calculate difference from linear shape factor distribution
dHb = Hb_bl - interp1(s_bl([1 end]),Hb_bl([1 end]),s_bl);

% Calculate integrated areas
q_1 = s_bl < s_bl_mid; q_2 = s_bl >= s_bl_mid;
p.aHb_1 = trapz(s_bl(q_1),dHb(q_1)); p.aHb_2 = trapz(s_bl(q_2),dHb(q_2));

% Record trailing edge shape factor value
p.Hb_te = Hb_bl(end);

% Record loss
p.loss = Polarx.omega(j);
p.loss_fr = Polarx.omega(j)-Polarx.omegv(j);

p.Mout = Polarx.mout;



%% Plot MISES results for the current section

% Check if plotting flag is set
if plot_stuff == 1 || plot_stuff == 2  
    % Open a figure window if requried or switch to one that is already open
    if exist('h','var') == 0 || isempty(h) == 1 
        
        % Open and resize figure window
        h.window = figure(); set(h.window,'Position',[1 12 960 800])
% % 		set(h.window,'Position',[1 62 1280 895])
%         
% %         % Open subplot windows
% h.ax(1) = axes;
        h.ax(1) = axes('position',[0.08 0.62 0.24 0.32]); 
        h.ax(2) = axes('position',[0.38 0.62 0.24 0.32]);
        h.ax(3) = axes('position',[0.68 0.62 0.24 0.32]); 
        h.ax(4) = axes('position',[0.08 0.22 0.24 0.32]);
        h.ax(5) = axes('position',[0.38 0.22 0.24 0.32]); 
        h.ax(6) = axes('position',[0.68 0.22 0.24 0.32]);
        h.ax(7) = axes('position',[0.08 0.04 0.39 0.11],'visible','off'); 
        h.ax(8) = axes('position',[0.53 0.04 0.39 0.11],'visible','off'); 
        
    elseif exist('h','var') ~= 0
        figure(h.window);
    else
        h = [];
    end
    
    % Turn on grid and boxes
    for n = 1:length(h.ax); axes(h.ax(n)); hold on; grid on; box on; end;
    
    % Switch to pressure distribution plot
%     axes(h.ax(1)); title('Pressure Distribution'); xlabel('Chord'); ylabel('Pressure Coefficient');
	axes(h.ax(1)); 
	title('Isentropic Mach number');  xlabel('Chord'); ylabel('Isentropic Mach'); ylim([0,1.4]);
    
    % Plot pressure distribution with peak suction
%     plot(s_1,-Cp_1,'-', 'Color', col); 
%     plot(p.s_Cp_max,interp1(s,-Cp,p.s_Cp_max,'pchip'),'.', 'Color', col); 
%     plot(s_2,-Cp_2,'-', 'Color', col);
	if exist('g','var')==1
		[blade]=ts_plot_2D(g,'M_is',nj_mid,[0,1.4],0,'Linestyle','--','Color',col);
		f=gcf;
		close(f);
		s1 = [0 ; cumsum(sum(diff(blade.xrt_1,1,1).^2,2).^0.5)]; s1 = s1 / max(s1);
		s2 = [0 ; cumsum(sum(diff(blade.xrt_2,1,1).^2,2).^0.5)]; s2 = s2 / max(s2);
		axes(h.ax(1))
% 		plot(s1,blade.var_1,'Color',col,'Linestyle','--');
% 		plot(s2,blade.var_2,'Color',col,'Linestyle','--');
	end
	axes(h.ax(1))
    plot(s_1,Mis_1,'-', 'Color', col); 
	p.s_1 = s_1;
    p.Mis_1 = Mis_1;
%     plot(p.s_Cp_max,interp1(s,-Cp,p.s_Cp_max,'pchip'),'.', 'Color', col); 
    plot(s_2,Mis_2,'-', 'Color', col);
    p.s_2 = s_2;
    p.Mis_2 = Mis_2;

%     Resize plot
    axis auto; v = axis; axis([0 1 0 1.4]);
    
    % Switch to shape factor plot
    axes(h.ax(2)); axis([0 1 1.5 3.5]); title('Shape Factor Distribution');
    xlabel('Chord'); ylabel('Boundary Layer Shape Factor');
    
    % Plot shape factor with ideal linear distribution
    plot(s_1,Hb_1,'-', 'Color', col);
% 	plot([s_bl_min s_bl_max],interp1(s_1,Hb_1,[s_bl_min s_bl_max]),'--', 'Color', col);
% 	bl=ts_bl_prop(g,1,50,0,1);
% 	plot(bl.s(:,nj_mid)/max(bl.s(:,nj_mid)),bl.H(:,nj_mid),'--b')
%     plot([s_bl_min s_bl_max],interp1(s_1,Hb_1,[s_bl_min s_bl_max]),'--', 'Color', col);
    
    %% Switch to leading edge plot
    axes(h.ax(3)); axis equal; axis(0.05*c.tchord * [-0.1 0.22 -0.1 0.2]); 
    title('Incidence'); xlabel('Axial Coordinate'); ylabel('Tangential Coordinate');
	
	% Plot stagnation point on the leading edge and the stagnation streamline
	plot(Polarx.xb{1} - c.mt_le(1),Polarx.yb{1} - c.mt_le(2),'-', 'Color', col,'Linewidth',1.0)
	plot(0,0,'.', 'Color', col)
	
	plot([cosd(phi_stag) ; 0 ] + mt_stag(1)-c.mt_le(1), ...
		[sind(phi_stag) ; 0 ] + mt_stag(2)-c.mt_le(2),'-.', 'Color', col)
	plot([c.mt_le_cen(1)- c.mt_le(1) ; mt_stag(1)-c.mt_le(1)], ...
		[c.mt_le_cen(2)- c.mt_le(2) ; mt_stag(2)-c.mt_le(2) ],'-.', 'Color', col,'Linewidth',1.0)
	plot([0 c.mt_le_cen(1)- c.mt_le(1)], [0 c.mt_le_cen(2)- c.mt_le(2)],'--', 'Color', col,'Linewidth',1.0)
	
	mt_stag_stream = [Idat.x(1:Idat.ninl(2),end) Idat.y(1:Idat.ninl(2),end)];
	mt_stag = [Idat.x(Idat.ninl(2),end) Idat.y(Idat.ninl(2),end)];
	
	plot(mt_stag_stream(:,1)- c.mt_le(1),mt_stag_stream(:,2)- c.mt_le(2),'-','Color',col,'Linewidth',1.0);
	plot(mt_stag(1)- c.mt_le(1),mt_stag(2)- c.mt_le(2),'x','Color',col,'Linewidth',1.0)
	plot(c.mt_le_cen(1)- c.mt_le(1),c.mt_le_cen(2)- c.mt_le(2),'o','Color',col,'Linewidth',1.0);
	

%     % Plot stagnation point on the leading edge and the stagnation streamline
%     plot(Polarx.xb{1} - c.mt_le(1),Polarx.yb{1} - c.mt_le(2),'-', 'Color', col)
%     plot(0,0,'.', 'Color', col)
%     plot([cosd(psi_stag + 90) ; 0 ; -1*cosd(psi_stag + 90)] + mt_stag(1)-c.mt_le(1), ...
%         [sind(psi_stag + 90) ; 0 ; -1*sind(psi_stag + 90)] + mt_stag(2)-c.mt_le(2),'-', 'Color', col)
%     plot([0 cosd(c.chi_le)], [0 sind(c.chi_le)],'--', 'Color', col)

	%%
	
    % Switch to momentum thickness plot or loss loop based on Polarx length
    axes(h.ax(4));
    
    if length(Polarx.binl) == 1
    title('Momentum Distribution'); xlabel('Chord'); ylabel('Momentum Thickness');
% 	title('Suction surface angle'); xlabel('Chord'); ylabel('Suction surface angle');
% 	plot(scl_ss,psi_ss,'-','Color',col)
    
%     Plot boundary layer development
    plot(s_1,Th_1.*Hb_1*c.r_le/job.tchord,'-','Color',col)
%     plot(s_2,Th_2.*Hb_2*c.r_le,'--','Color',col)

    else

    title('Loss Loop'); xlabel('Incidence relative to design flow angle'); ylabel('Profile loss');
    plot(Polarx.binl-Ises.binl,Polarx.omega,'Color',col)
    end
    
    % Resize plot
    axis auto;

% 	bl=ts_bl_prop(g,1,50,0,0);
% 	plot(bl.s(:,nj_mid)/max(bl.s(:,nj_mid)),bl.th(:,nj_mid),'-r')
    
    % Plot transition locations
    xle = min(Polarx.xb{1}); xte = (Polarx.xb{1}(1) + Polarx.xb{1}(end))*0.5;
    strans = (Polarx.xtr(1,j) - xle) / (xte - xle);
%     plot(strans,interp1(s_1,Th_1,strans),'.','Color',col)
    strans = (Polarx.xtr(2,j) - xle) / (xte - xle);
%     plot(strans,interp1(s_2,Th_2,strans),'.','Color',col)
    
    
    % Plot non-dimensional camber line
    axes(h.ax(5)); 
%     axis([0 1 0 1]); 
    title('Camber Line'); xlabel('Chord'); ylabel('Camber');
    cam_spline = 1;
    cam_param = bl_construct_camber(job,[],cam_spline);
    plot(cam_param.s,1-cam_param.cam,'-','Color',col);
    plot(cam_param.s,cam_param.chi(1)-cam_param.chi,'--','Color',col);
    xlim([0 1])

    
    % Plot the whole blade geometry
    axes(h.ax(6)); axis equal;
    title('Geometry'); xlabel('Axial Coordinate'); ylabel('Tangential Coordinate');
%     plot(c.xrt_cam(:,1),c.xrt_cam(:,2),'-','Color',col)
%     plot(c.xrt_cl(:,1),c.xrt_cl(:,2),'-','Color',col)
	if exist('job','var')==1
	plot(xrt(:,1),xrt(:,2),'--','Color',col)
	plot((Polarx.xb{1})*c.r_le,(Polarx.yb{1})*c.r_le,'-', 'Color', col,'Linewidth',1.0)
	end
    
    % Find how many lines there are on the camber plot
    np = length(findobj(h.ax(5),'Type','line'));
    
    % Text names
    varnames_1 = {'psi_stag' 'M_max'  's_Cp_max' 'Alpha' 'loss'};
    texnames_1 = {'\Psi (local incidence ^{o})' 'M_{peak}' 's_M_{peak}' '\alpha_{out} (^{o})' '\omega (stagnation pressure loss)'};
%     varnames_2 = {'s_thick_max' 'chi_te' 'chi_le' 's_c'};
%     varnames_2 = {'oscosa' 'chi_te' 'chi_le' 's_c'};
    varnames_2 = {'AtA1' 'chi_te' 'chi_le' 's_c' 't_c'};
    texnames_2 = {'A_{throat}/A_{inlet}' '\chi_{te} (^{o})' '\chi_{le} (^{o})' 's/c' 't_{max}/c'};    
     
    % Print text titles
    if np == 1
        
        % Print flow names
        axes(h.ax(7));
        for v = 1:length(texnames_1); 
            text(0,v,[texnames_1{v} ' = '],'FontSize',11,'Color',[0 0 0])
        end
        
        % Print geometry names
        axes(h.ax(8));
        for v = 1:length(texnames_2); 
            text(0,v,[texnames_2{v} ' = '],'FontSize',11,'Color',[0 0 0])
        end        
    end
    
    % Print flow values
    axes(h.ax(7)); axis([0 np+1 1 length(texnames_1)]);
    for v = 1:length(varnames_1); 
        text(np,v,sprintf('%+.3f',p.(varnames_1{v})),'FontSize',11,'Color',col)
    end
    
    % Print all geometry values
    axes(h.ax(8)); axis([0 np+1 1 length(texnames_2)]);
    for v = 1:length(varnames_2); 
        text(np,v,sprintf('%+.3f',c.(varnames_2{v})),'FontSize',11,'Color',col)
    end    
    
    %% plot contour lines
    if plot_stuff == 2
    H=figure; grid on; hold on;
    H.Position = [471.0000 143 331.0000 603.7772];
    % determine contour levels
    M_contours = [0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.05 1.1 1.15 1.2 1.25 1.30 1.35 1.40];
    [C_pitch,f]=contour(Idat.x_cell.*c.r_le,Idat.y_cell.*c.r_le+pitch,Idat.M,M_contours,'Color',[0.25,0.25,0.25]);
    p.C_pitch = C_pitch;
    % plot sonic contour in red
    [C,f]=contour(Idat.x_cell.*c.r_le,Idat.y_cell.*c.r_le+pitch,Idat.M,[1.0 1.0],'Color',[1,0,0],'Linewidth',1.1);
    [C,f]=contour(Idat.x_cell.*c.r_le,Idat.y_cell.*c.r_le,Idat.M,[1.0 1.0],'Color',[1,0,0],'Linewidth',1.1);
    % plot rest of contours
    [C,f]=contour(Idat.x_cell.*c.r_le,Idat.y_cell.*c.r_le,Idat.M,M_contours,'Color',[0.25,0.25,0.25]);
    p.C = C;
  
    % make plot limits
    c_x = max(c.xrt(:,1)) - min(c.xrt(:,1));
    axis equal; 
    xlim([min(c.xrt(:,1))-0.5*c_x,max(c.xrt(:,1))+0.3*c_x])
    ylim([min(c.xrt(:,2))-1.30*pitch,max(c.xrt(:,2))+0.30*pitch])
    % manually labe or automatically label
    %     clabel(C,h,'LabelSpacing',100,'FontSize',9);
    clabel(C,f,'LabelSpacing',400,'FontSize',9);
%     clabel(C,f,'manual','FontSize',9);
    % plot blade profile
    plot(c.xrt(:,1),c.xrt(:,2),'-','Color',[0,0,0],'Linewidth',1.3)
    plot(c.xrt(:,1),c.xrt(:,2)-pitch,'-','Color',[0,0,0],'Linewidth',1.3)
    plot(c.xrt(:,1),c.xrt(:,2)+pitch,'-','Color',[0,0,0],'Linewidth',1.3)
    % plot boundary layer
    plot(bl_stream(:,1),bl_stream(:,2),'Color',[0,0,1],'Linestyle','--','Linewidth',1.3)
    plot(bl_stream(:,1),bl_stream(:,2)-pitch,'Color',[0,0,1],'Linestyle','--','Linewidth',1.3)
    % plot throat location
    plot(xrt_throat(:,1).*pitch,xrt_throat(:,2).*pitch,'Color',[0,0.5,0],'Linestyle','-.','Linewidth',1.3)
    plot(xrt_throat(:,1).*pitch,xrt_throat(:,2).*pitch-pitch,'Color',[0,0.5,0],'Linestyle','-.','Linewidth',1.3)
    plot(xrt_throat_bl(:,1).*pitch,xrt_throat_bl(:,2).*pitch,'Color',[0,0.5,0],'Linestyle',':','Linewidth',1.3)
    plot(xrt_throat_bl(:,1).*pitch,xrt_throat_bl(:,2).*pitch-pitch,'Color',[0,0.5,0],'Linestyle',':','Linewidth',1.3)
    % create label close to throat
    text(double(0.70*xrt_throat(1,1)+0.30*xrt_throat(end,1)).*pitch+0.025*c_x,double(0.30*xrt_throat(end,2)+0.70*xrt_throat(1,2)).*pitch-0.95*pitch,['A_{t}/A_{1}=' num2str(round(100*p.AtA1)/100)],'Color',[0,0.5,0],'FontName','Arial');
    text(double(0.20*xrt_throat_bl(1,1)+0.80.*xrt_throat_bl(end,1)).*pitch+0.025*c_x,double(0.80*xrt_throat_bl(end,2)+0.20*xrt_throat_bl(1,2)).*pitch-0.95*pitch,['A_{t}/A_{1}=' num2str(round(100*p.AtA1_bl)/100)],'Color',[0,0.5,0],'FontName','Arial');
    plot([(0.50*xrt_throat_bl(1,1)+0.50.*xrt_throat_bl(end,1)).*pitch,(0.20*xrt_throat_bl(1,1)+0.80.*xrt_throat_bl(end,1)).*pitch+0.15*c_x],[(0.50*xrt_throat_bl(1,2)+0.50.*xrt_throat_bl(end,2)).*pitch-pitch,(0.20*xrt_throat_bl(1,2)+0.80.*xrt_throat_bl(end,2)).*pitch+0.025*c_x-pitch],'Color',[0,0,0])
     plot([(0.55*xrt_throat(1,1)+0.45.*xrt_throat(end,1)).*pitch,(0.70*xrt_throat(1,1)+0.30.*xrt_throat(end,1)).*pitch+0.025*c_x],[(0.55*xrt_throat(1,2)+0.45.*xrt_throat(end,2)).*pitch-pitch,(0.70*xrt_throat(1,2)+0.30.*xrt_throat(end,2)).*pitch+0.025*c_x-0.95*pitch],'Color',[0,0,0])
    tightfig;
    axis off;
    
    end

% elseif plot_stuff == 2;
% 	 if exist('h','var') == 0 || isempty(h) == 1
% 		 h.window=figure(); hold on; grid on; axis on;
% 	 end
% 	figure(h.window);	
% 	h.ax(1)=subplot(1,2,1); grid on; hold on;
% 	h.ax(2)=subplot(1,2,2); grid on; hold on;
% %% plot isentropic Mach number
% 	axes(h.ax(2)); 
% %  	title('Isentropic Mach number');  
% 	xlabel('s/s_{o}'); ylabel('Isentropic Mach number'); ylim([0,1.4]); axis square;
% % 	xlabel('Non-dimensional chord'); ylabel('Pressure coefficient (-C_{p})');
% 	hold on;
%     plot([s_1;flip(s_2)],[Mis_1;flip(Mis_2)],'-', 'Color', col,'Linewidth',1); 
% % 	
% %% 	% Plot shape factor with ideal linear distribution
% % 	plot([s_1;flip(s_2)],[-Cp_1;flip(-Cp_2)],'-', 'Color', col); 
% % 	axes(h.ax(2));
% % 	xlabel('Non-dimensional chord'); ylabel('Boundary Layer Shape Factor'); axis square;
% % 
% % 	plot(s_1,Hb_1,'-', 'Color', col);
% % 	% 	plot([s_bl_min s_bl_max],interp1(s_1,Hb_1,[s_bl_min s_bl_max]),'--', 'Color', col);
% % 	% 	bl=ts_bl_prop(g,1,50,0,1);
% % 	% 	plot(bl.s(:,nj_mid)/max(bl.s(:,nj_mid)),bl.H(:,nj_mid),'--b')
% % 	plot([s_bl_min s_bl_max],interp1(s_1,Hb_1,[s_bl_min s_bl_max]),'--', 'Color', col,'Linewidth',1);
% % 	ylim([1,3])
% %% plot suction surface angle
% 	axes(h.ax(1)); xlabel('s/s_{o}'); ylabel('Suction surface angle');
% 	plot(scl_ss,psi_ss,'-','Color',col)
% % 	plot(scl_ps,psi_ps,'--','Color',col)
% 	axis square
% %% Plot non-dimensional camber line
% %     axes(h.ax(2));  xlabel('x/c_{x}'); ylabel('\kappa/\kappa_{tot}'); axis square
% % % 	axis([0 1 0 1]);
% % %     plot(c.s_cl,c.cam,'-','Color',col);
% % 	if exist('job','var')==1
% % 		scl = [0 ; cumsum(sum(diff(mca.xrt_ss,1,1).^2,2).^0.5,1)];
% % 		size(mca.s_cl)
% % 		size(mca.cam)
% % 	plot(scl/max(scl),mca.cam,'Color',col,'Linestyle','-');
% % 	end
% % % 	axis square
% % % 	
% %% Plot local incidence
% %     % Switch to leading edge plot
% %     axes(h.ax(2)); axis equal; axis(0.5*c.tchord * [-0.1 0.22 -0.1 0.2]);
% %     xlabel('Axial Coordinate'); ylabel('Tangential Coordinate');
% %     
% % %     % Plot stagnation point on the leading edge and the stagnation streamline
% %     plot(Polarx.xb{1} - c.mt_le(1),Polarx.yb{1} - c.mt_le(2),'-', 'Color', col,'Linewidth',1.0)
% % %     plot(0,0,'.', 'Color', col)
% % % 	
% % %     plot([cosd(phi_stag) ; 0 ] + mt_stag(1)-c.mt_le(1), ...
% % %         [sind(phi_stag) ; 0 ] + mt_stag(2)-c.mt_le(2),'-.', 'Color', col)
% % 	plot([c.mt_le_cen(1)- c.mt_le(1) ; mt_stag(1)-c.mt_le(1)], ...
% % 		[c.mt_le_cen(2)- c.mt_le(2) ; mt_stag(2)-c.mt_le(2) ],'-.', 'Color', col,'Linewidth',1.0)
% %     plot([0 c.mt_le_cen(1)- c.mt_le(1)], [0 c.mt_le_cen(2)- c.mt_le(2)],'--', 'Color', col,'Linewidth',1.0)
% % % 	
% % 	mt_stag_stream = [Idat.x(1:Idat.ninl(2),end) Idat.y(1:Idat.ninl(2),end)];
% % 	mt_stag = [Idat.x(Idat.ninl(2),end) Idat.y(Idat.ninl(2),end)];
% % 
% % 	plot(mt_stag_stream(:,1)- c.mt_le(1),mt_stag_stream(:,2)- c.mt_le(2),'-','Color',col,'Linewidth',1.0);
% % 	plot(mt_stag(1)- c.mt_le(1),mt_stag(2)- c.mt_le(2),'x','Color',col,'Linewidth',1.0)
% % 	plot(c.mt_le_cen(1)- c.mt_le(1),c.mt_le_cen(2)- c.mt_le(2),'o','Color',col,'Linewidth',1.0);
% %% plot local incidence against total camber
% % 	if exist('job','var')==1
% % 		 axes(h.ax(2));
% % 		Polarx.binl(j)
% % 		c.chi_le
% % 		job.tot_cam
% % 	plot(Polarx.binl(j)-c.chi_le,job.tot_cam,'x','Color',col);
% % 	end
% % 	xlabel('\alpha_{1}-\chi_{1}'); ylabel('Total camber \chi_{2}-\chi_{1}');

elseif plot_stuff == 3;
	if exist('h','var') == 0 || isempty(h) == 1
		h=figure(); hold on; grid on; axis on;
	end
	figure(h);
		xlabel('Non-dimensional chord'); ylabel('Isentropic Mach number'); ylim([0,1.4]); 
%         axis square;
% 	xlabel('Non-dimensional chord'); ylabel('Pressure coefficient (-C_{p})');
	hold on;
    plot([s_1;flip(s_2)],[Mis_1;flip(Mis_2)],'-', 'Color', col,'Linewidth',1.5);
end    

end