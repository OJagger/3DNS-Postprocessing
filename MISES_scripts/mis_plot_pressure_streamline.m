function [h] = mis_plot_pressure_streamline(directory,h,col)
% function to plot isentropic pressure streamline and suction surface streamline like in D.Lefas Paper
%
%   [p,h] = MIS_PLOT_SECTION(directory,h,col,plot_stuff,t)
%
%   directory - string of output file directory
%   h - optional figure handle
%   col - optional colour RGB vectors
%   plot_stuff - 0 or 1 fgeditor showing working
%   t - optional blade section information
%   p - output data structure

% % Default to show plot
% if exist('plot_stuff','var') == 0
%     plot_stuff = 1;
% end

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

% Find the closest point to peak suction
s = Polarx.s{j}(:,1) / Polarx.s{j}(Polarx.iteb(1),1);


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
    return
end


if exist('h','var') == 0 || isempty(h) == 1 
    h.pressure = figure(); grid on; hold on; 
    h.dPo = figure(); grid on; hold on; 
end

bl_inc = 1;
[rad_contr,AtA1,AtA1_rhoV,o_s,rad_contr_o,xrt_throat] = mis_radial_contraction(directory,0,bl_inc);
display(['At/A1 = ' num2str(AtA1)])

%% find isnetropic streamline
figure(h.pressure);
j_is=max(find(Idat.dPo_Po(end,:)>0,2))
if isempty(j_is) == 1
   j_is = size(Idat.dPo_Po,2);
end
% find center of loss core
[~,j_loss]=min(Idat.dPo_Po(end,:));
j_mid = round(size(Idat.dPo_Po,2)./2);

% find inlet/le and outlet/te plane to determine the amount s the streamline has
% travelled
load([directory 'section.mat']);
[~,i_le]=min(abs(Idat.x(:,1)));
[~,i_te] = min(abs(Idat.x(:,end)-c.m_chord));
Idat.s_le=squeeze(Idat.s(i_le,1));
Idat.s_te=squeeze(Idat.s(i_te,1));


% plot((Idat.s(:,j_is)-Idat.s_le(j_is))./((Idat.s_te(j_is)-Idat.s_le(j_is))),Idat.P_Po1(:,j_is),'--','Color',col)
% plot((Idat.s(:,j_mid)-Idat.s_le(j_mid))./((Idat.s_te(j_mid)-Idat.s_le(j_mid))),Idat.P_Po1(:,j_mid),'.-','Color',col)
% plot((Idat.s(:,j_loss)-Idat.s_le(j_loss))./((Idat.s_te(j_loss)-Idat.s_le(j_loss))),Idat.P_Po1(:,j_loss),'-','Color',col)

plot((Idat.s(:,j_is)-Idat.s_le)./(Idat.s_te-Idat.s_le),Idat.P_Po1(:,j_is),'--','Color',col)
plot((Idat.s(:,j_mid)-Idat.s_le)./(Idat.s_te-Idat.s_le),Idat.P_Po1(:,j_mid),'.-','Color',col)
plot((Idat.s(:,j_loss)-Idat.s_le)./(Idat.s_te-Idat.s_le),Idat.P_Po1(:,j_loss),'-','Color',col)


% plot(Idat.x(i_le(j_is),j_is),Idat.P_Po1(i_le(j_is),j_is),'x','Color',col)
% plot(Idat.x(i_te(j_is),j_is),Idat.P_Po1(i_te(j_is),j_is),'x','Color',col)

legend('Isentropic','Mid-streamline','Loss core')
xlabel('Streamwise distance')
ylabel('P/Po')
% xlim([0 , 1])

%% plot loss core
figure(h.dPo);
plot(Idat.dPo_Po(end,:),(Idat.y_cell(end,:)-min(Idat.y_cell(end,:)))./(max(Idat.y_cell(end,:))-min(Idat.y_cell(end,:))),'Color',col)
p1=plot(Idat.dPo_Po(end,j_is),(Idat.y_cell(end,j_is)-min(Idat.y_cell(end,:)))./(max(Idat.y_cell(end,:))-min(Idat.y_cell(end,:))),'x','Color',col);
p2=plot(Idat.dPo_Po(end,j_mid),(Idat.y_cell(end,j_mid)-min(Idat.y_cell(end,:)))./(max(Idat.y_cell(end,:))-min(Idat.y_cell(end,:))),'sq','Color',col);
p3=plot(Idat.dPo_Po(end,j_loss),(Idat.y_cell(end,j_loss)-min(Idat.y_cell(end,:)))./(max(Idat.y_cell(end,:))-min(Idat.y_cell(end,:))),'o','Color',col);
legend([p1,p2,p3],'Isentropic','Mid-streamline','Loss core')
xlabel('\Delta Po/Po')
end