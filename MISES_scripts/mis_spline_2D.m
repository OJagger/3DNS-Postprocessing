% pass a spline through the suction surface Mach number distribution
%Written by Demetrios Lefas (dl467), 2017
function [Area1,Area2,M_le,M_max,x_max,x_max_orig] = mis_spline_2D(directory,M_factor,plot_stuff,xrange)

directory = strrep(directory,'TURBOSTREAM','MISES');
load([directory 'section.mat'])

%% extract parameters from polarx.mises file

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
[~,j] = min(abs(Polarx.binl-Ises.binl));

% Record the outlet flow angle
p.Alpha = atand(Polarx.sout(j));

% Find the closest point to peak suction
Cp = Polarx.cp{j}(:,1); s = Polarx.s{j}(:,1) / Polarx.s{j}(Polarx.iteb(1),1);
Mis = Polarx.mn{j}(:,1);
q = s > 0.01 & s < 0.99; Cp = Cp(q); s = s(q); Mis = Mis(q);
[~, i] = min(Cp); i = max(i,3);

% Fit a polynomial through points near peak suction to find the maximum more accurately
s_temp = linspace(s(i-2),s(i+2),1000);
Q = polyval(polyfit(s(i-2:i+2),Cp(i-2:i+2),3),s_temp);
[~,i_max] = min(Q); p.s_Cp_max= s_temp(i_max);
M = polyval(polyfit(s(i-2:i+2),Mis(i-2:i+2),3),s_temp);
[~,i_max] = min(M); p.s_M_max= s_temp(i_max);

% Read in the grid coodinates
if exist([directory 'idat.mises_01'],'file') ~= 0
    Idat = mis_read_idat('mises_01',directory);
else
    Idat = mis_read_idat('mises',directory);
end

% Calculate the relative angle from the leading edge centre
mt_stag = [Idat.x(Idat.ninl(2),end) Idat.y(Idat.ninl(2),end)];
phi_stag = atan2(c.mt_le_cen(2) - mt_stag(2),c.mt_le_cen(1) - mt_stag(1)) * 360 / (2*pi);
mt = c.mt(c.mt(:,1) < c.mt_le_cen(1),:);
phi = atan2(c.mt_le_cen(2) - mt(:,2),c.mt_le_cen(1) - mt(:,1)) * 360 / (2*pi);

% Calculate the angle the stagnation point makes with the leading edge
[phi,i] = unique(phi); mt = mt(i,:);
dtdx = grad_mg(mt(:,1),mt(:,2));
psi = atand(dtdx);

% Correct the angle by the leading edge metal angle
p.psi_stag = 90 - (c.chi_le - interp1(phi,psi,phi_stag,'pchip'));
psi_stag = interp1(phi,psi,phi_stag,'pchip');

% Record the boundary layer parameters
q_1 = Polarx.ileb(1):Polarx.iteb(1);
q_2 = Polarx.ileb(2):Polarx.iteb(2);
s_1 = Polarx.s{j}(q_1,1); s_1 = s_1 / max(s_1);
s_2 = Polarx.s{j}(q_2,2); s_2 = s_2 / max(s_2);

Cp_1 = Polarx.cp{j}(q_1,1); Cp_2 = Polarx.cp{j}(q_2,2);
Mis_1 = Polarx.mn{j}(q_1,1); Mis_2 = Polarx.mn{j}(q_2,2);
Hb_1 = Polarx.hbar{j}(q_1,1);
Th_1 = Polarx.th{j}(q_1,1); Th_2 = Polarx.th{j}(q_2,2); 

% Record the shape factor profile parameters
p.s = s_1; p.Hb = Hb_1;

%%

x_ss = s_1;
var_ss = Mis_1;


% find maximum (remove leading edge)
var_ss_max = var_ss(x_ss>0.05);
x_ss_max = x_ss(x_ss>0.05);
[var_max,i_max]=max(var_ss_max);
x_max=x_ss_max(i_max);
x_max_orig = x_max;

% find M_max 
[M_max,~]=max(var_ss_max);

if M_max<1.05
if exist('xrange','var')==0
if x_max>0.47 % otherwise knot sequence isn't increasing
    x_max=0.45;
end
if x_max<0.08
    x_max = 0.15; % otherwie it creates a too sudden expansion
end
else
    if x_max>xrange(2)
        x_max = xrange(2) - 0.05;
	end
    if x_max<xrange(1);
        x_max = xrange(1) + 0.05;
	end
end
end

knots = [x_max, 2*x_max, 2*x_max]; order = 4;
% Add on starting and end knots
knots = [0.01*ones(1,order)  knots 0.98*ones(1,order)];

d = order - 1; %degree
t_star=zeros(1,length(knots)-(d+1));
for j=1:length(knots)-(d+1)
    for m=1:d
    t_star(j) = t_star(j) + knots(j+m);
    end
end
t_star=t_star/d;

var_ss = var_ss(x_ss>0.005);
x_ss = x_ss(x_ss>0.005);


% in case it goes round the leading edge
[~,I]=min(x_ss);
x_ss = x_ss(I:end);
var_ss = var_ss(I:end);

% in case it goes round the trailing edge
[~,J]=max(x_ss);
x_ss = x_ss(1:J);
var_ss = var_ss(1:J);
% keep only unique values
[x_ss,K] = unique(x_ss);
var_ss = var_ss(K);

coef = zeros(1,length(t_star));


for i= 1:length(coef)
coef(i) = interp1(x_ss,var_ss,t_star(i),'linear');
end



%% M_le is the inlet Mach number
M_le = f.M;
% coef(1) = M_le;
%%

%% NEED TO PUT JOB.PSI_TARGET in at some point
psi_target =0;
%%
M_le
% % in case of a huge leading edge spike
if coef(1) > 1.2*M_le || coef(1) > 1 || abs(p.psi_stag-psi_stag)>4 || coef(1)<M_le  % so that it can combat huge spikes
%% I have changed this
%         coef(1)=0.75*M_le+0.25*coef(1);
% elseif coef_ps_le > M_le*1.1
%    coef(1)=0.75*M_le+0.25*coef_ps_le;
    coef(1) = M_le+0.08;

end

% coef(3)=1.05*coef(3); % so that spline approximation reaches closer to the peak
if isempty(M_factor) == 0
coef(3)=1.03*M_factor*M_le;
else
    if M_max > 1.30
        M_max = 1.30;
	end
    coef(3) = 1.03*M_max;    
end
% coef(2)=(coef(1)+coef(3))/2;
coef(2)=(0.75*coef(2)+0.25*coef(3));
if coef(4)>coef(3) % in case the maximum is now coef(4)
coef(4)=(0.75*coef(3)+0.25*coef(5));
else
   coef(4)=(0.15*coef(3)+0.85*coef(4)); 
end


%% split into two splines (pre and after shock)
if M_max>1.05 
%     && x_max == x_max_orig

    I_max=0;
    % determine whether is has too bumps (reacceleration - have a lot of camber at the front)
	[~,i_min]=min(abs(var_ss_max-1));
    x_ss_M_1=x_ss_max(i_min);
	if x_max_orig < x_ss_M_1-0.15
		x_max_orig = x_ss_M_1-0.15;
		I_max =1;
	end
	
    alpha=0.675;
    beta=0.34;
    lamda3=(3*alpha-2)+3*(1-alpha)*x_max_orig;
    lamda2=(3*beta-1)+3*(1-beta)*x_max_orig-lamda3;
    if x_max_orig > 1 - 1/(3*beta)
    lamda1=3*x_max_orig-lamda2-lamda3;
    else
        lamda1=0.01;
%         lamda2=0;
%         lamda3=0;
    end
    
    if lamda1 < 0.02 % if coeficient is too close too leading edge just remove because it doesn't contribute to the solution
    lamda1 = 0.01;
    end
    
    
    % cubic B-spline
    knots = [lamda1 lamda2 lamda3];
    knots = [0.01*ones(1,order)  knots 0.98*ones(1,order)];
    
    d = order - 1; %degree
    t_star=zeros(1,length(knots)-(d+1));
    for j=1:length(knots)-(d+1)
        for m=1:d
            t_star(j) = t_star(j) + knots(j+m);
        end
    end
    t_star=t_star/d;

    coef = zeros(1,length(t_star));

    for i= 1:length(coef)
        coef(i) = interp1(x_ss,var_ss,t_star(i),'linear');
    end
    
    
    if isempty(M_factor) == 0
        if coef(4)<coef(3)
        coef(3)=M_factor*M_le;
        else % in case the maximum is now coef(4)
           coef(4) = M_factor*M_le;
           coef(3)=(0.25*coef(4)+0.75*coef(3));
        end
    else
        if coef(4)<coef(3) && I_max == 0
        coef(3) = 0.98*M_max;
        else % in case the maximum is now coef(4)
           coef(4) = 0.98*M_max;
           coef(3)=(0.25*coef(4)+0.75*coef(3));
        end
    end
    % coef(2)=(coef(1)+coef(3))/2;
    coef(2)=(0.75*coef(2)+0.25*coef(3));

    % % in case of a huge leading edge spike
%     if coef(1) > 1.2*M_le  % so that it can combat huge spikes
%         coef(1)=0.25*M_le+0.75*coef(1);
%     elseif coef_ps_le > M_le*1.1
%         coef(1)=0.75*M_le+0.25*coef_ps_le;
%     end
if coef(1) > 1.2*M_le || coef(1) > 1 || abs(p.psi_stag-psi_target)>4 || coef(1)<M_le  % so that it can combat huge spikes
    %% I have changed this
    %         coef(1)=0.75*M_le+0.25*coef(1);
    % elseif coef_ps_le > M_le*1.1
    %    coef(1)=0.75*M_le+0.25*coef_ps_le;
    
    coef(1) = M_le*1.04;
    
end

    % drop subsonic value right after shock (too make it more straight)
    coef(5) = 0.90*coef(5);
end


knots;
t_star;
coef;

spline_var = spmak(knots,coef);

% calculate spline values within range
var_ss = var_ss(x_ss>0.01 & x_ss<0.98);
x_ss = x_ss(x_ss>0.01 & x_ss<0.98);

var_spline = fnval(spline_var,x_ss);

delta_var = abs(var_spline - var_ss);
for i=1:length(delta_var)-1
   delta_var_avg(i) = (delta_var(i)+delta_var(i+1))/2;
   diff(i) = x_ss(i)-x_ss(i+1);
   x_avg(i) = (x_ss(i)+x_ss(i+1))/2;
end

% delta_var_avg
% x_avg

Area1 = abs(sum(delta_var_avg(x_avg<t_star(3)).*diff(x_avg<t_star(3))));
Area2 = abs(sum(delta_var_avg(x_avg>t_star(3)).*diff(x_avg>t_star(3))));
%%
% Area2 = 0;
%%
Area = sum(delta_var_avg.*diff);

if plot_stuff == 1
    figure;
    hold on
    plot(x_ss,var_ss)
    plot(x_ss,var_spline,'-r')
    plot(t_star,coef,'--x')
    
    grid on
    xlabel('$s/s_{max}$','Interpreter','latex')
    ylabel(['$M_{is}$'],'Interpreter','latex')
    xlim([0 1])
end