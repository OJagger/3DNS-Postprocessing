
% function that takes a mises datum blade and then tweaks the important
% variables that define Athroat/Ainlet (inlet blade metal angle, percentage camber, outlet blade metal angle, thickness, pitch-to-chord & Axt/Ax1)
function mis_sweep_datum(job,plot_stuff)

if exist('plot_stuff','var')==0
    plot_stuff = 0;
end

if plot_stuff == 1 
   h = []; % plot general data summary
   h1 = figure(); grid on; hold on;
   h2 = figure(); grid on; hold on;
   h3 = figure(); grid on; hold on;
end

job_datum = job;
plot_freeman = 0;

%% figure linestyle and color details details
color = [0,0,1;1,0,0;0,0.5,0;1,0,1;0.5,0,0.5;0.5,0.5,0.5;0,0,0;0.25,0.25,0.25;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0;0,0,0];
Linestyle = {'-','--',':','-.','-','-','-','-','-','-'};
marker = {'^','o','square','d','v','<','x','+','.'};

%% Extract file prefix
file_prefix = strrep(job_datum.rjm_directory,'/mnt/Disk1/dl467/Documents/Test_Cases/','');
file_prefix = strrep(file_prefix,'MISES_solutions/','');
file_prefix = strrep(file_prefix,'/MISES/','');
disp(['Datum blade is : ' file_prefix])


J = 1;
Legend = [];

job.new_datum_directory = ['/mnt/Disk1/dl467/Documents/Test_Cases/MISES_summary/' file_prefix '/'];
disp(['***** ' job.new_datum_directory ' *****'])
%% sweep through t/c
t_c_datum = job_datum.tchord/(2*pi*0.5*(job_datum.r_le +job_datum.r_te));
t_c_min = 0.01; t_c_max = 0.06;
dt_c = [-0.02 -0.01 0 0.01 0.02];
t_c_sweep = t_c_datum+dt_c;
t_c_sweep = t_c_sweep(t_c_sweep>t_c_min); t_c_sweep = t_c_sweep(t_c_sweep<t_c_max);

disp(['*** sweeping: t/c ****'])

for m = 1:length(t_c_sweep)
    h = [];
    % sweep through different perecentage cambers
    perc_cam_datum = job_datum.perc_cam;
    dperc_cam = [-0.2 -0.15 -0.1 -0.05 0 0.05 0.1 0.15 0.20];
    perc_min = 0; perc_max = 0.70;
    perc_cam = perc_cam_datum + dperc_cam;
    perc_cam = perc_cam(perc_cam>perc_min); perc_cam = perc_cam(perc_cam<perc_max);
    k = 1;
    for i = 1:length(perc_cam)
        
        job.perc_cam = perc_cam(i);
        directory = [job.new_datum_directory file_prefix '-' 'TC' num2str(round(100*t_c_sweep(m))) 'PC' num2str(round(100*perc_cam(i))) '/MISES/'];
        job.rjm_directory = directory;
        % Make directory if required
        if exist(directory,'dir') == 0
            mkdir(directory);
            disp([directory ': does not exist'])
        else
            [Polarx, Ises] = mis_read_polarx('mises',directory);
            if isstruct(Polarx) ~= 1
                disp([directory ' : solution did not converge'])
                continue
            else
                if plot_stuff == 1
                    % calculate Athroat
                    plot_blade = 0; ote = 0; dev = []; cam_spline = 1;
                    [xrt,xrt_cam,xrt_ps,xrt_ss]=bl_construct_section(job,plot_blade,ote,dev,cam_spline);
                    mca.xrt_ss = xrt_ss;
                    mca.xrt_ps = xrt_ps;
                    pitch = 2*pi*job.r_mid/job.N; mt_stag_stream = [];
                    [~,o_s(k),xrt_throat]=ts_throat_calc(mca.xrt_ss,mca.xrt_ps,pitch,mt_stag_stream/pitch,0);
                    oscosa(k)=o_s(k)/cosd(job.alpha_rel_inlet);
                    [rad_contr(k)] = mis_radial_contraction(directory);
                    AtA1(k) = oscosa(k).*rad_contr(k);
                    % plot mises data summary
                    [p,h] = mis_plot_section(directory,h,color(i,:),plot_stuff,job);
                    title([file_prefix '-' 'TCPC'])
                    % store Peak Mach number and loss (overall and freestream)
                    % Find the closest point to peak suction
                    M_max(k) = p.M_max_ss;
                    Yp(k) = p.loss;
                    Yp_fr(k) = p.loss_fr;
                end
                 if i ~= length(perc_cam)
                k = k+1;
                continue
                 end
            end
        end
        cd(directory)
        
        job.thick_max = t_c_sweep(m).*job_datum.tchord;

        % ADJUST FRACTION chord to betweeen covered passage and shock
        % (mainly for pitch-to-chord)
        [job]=ts_frac_chord_mises(job,[],2); % mode 2 - between shock and covered passage
        [q,h,job] = mis_run_section(directory,[],[0,0,1],0,1,0,job);
        if q.manual == 1
            disp(['****** Cannot run datum at : ' 'TC' num2str(round(100*t_c_sweep(m))) 'PC' num2str(round(100*perc_cam(i))) ' ******'])
        end
        save([directory 'job.mat'],'job')
        

    end
    
    Legend = [Legend ; {['TC = ' num2str(round(100*t_c_sweep(m))) '%']}];
    if plot_stuff == 1 && exist('AtA1','var')==1
        % PLOT FIGURES
        [h1,h2,h3]=plot_fig(h1,h2,h3,J,m,Legend);
        % clear all recalculated variables
        clear Yp Yp_fr M_max AtA1
        h = [];
    end

end



%% sweep through AxtAx1

J = J + 1;
k = 1;
job = job_datum;
AxtAx1_sweep = [0.95 0.975 1.025 1.05];
for m = 1:length(AxtAx1_sweep)
    h = [];
    % sweep through different perecentage cambers
    perc_cam_datum = job_datum.perc_cam;
    dperc_cam = [-0.2 -0.15 -0.1 -0.05 0 0.05 0.1 0.15 0.20];
    perc_min = 0; perc_max = 0.70;
    perc_cam = perc_cam_datum + dperc_cam;
    perc_cam = perc_cam(perc_cam>perc_min); perc_cam = perc_cam(perc_cam<perc_max);
    for i = 1:length(perc_cam)
        
        job.perc_cam = perc_cam(i);
        directory = [job.new_datum_directory file_prefix '-' 'Axt' num2str(round(100*AxtAx1_sweep(m))) 'PC' num2str(round(100*perc_cam(i))) '/MISES/'];
        job.rjm_directory = directory;
        % Make directory if required
        if exist(directory,'dir') == 0
            mkdir(directory);
            disp([directory ': does not exist'])
        else
            [Polarx, Ises] = mis_read_polarx('mises',directory);
            if isstruct(Polarx) ~= 1
                disp([directory ' : solution did not converge'])
                continue
            else
                if plot_stuff == 1
                    % calculate Athroat
                    plot_blade = 0; ote = 0; dev = []; cam_spline = 1;
%                     [xrt,xrt_cam,xrt_ps,xrt_ss]=bl_construct_section(job,plot_blade,ote,dev,cam_spline);
%                     mca.xrt_ss = xrt_ss;
%                     mca.xrt_ps = xrt_ps;
%                     pitch = 2*pi*job.r_mid/job.N; mt_stag_stream = [];
%                     [~,o_s(k),xrt_throat]=ts_throat_calc(mca.xrt_ss,mca.xrt_ps,pitch,mt_stag_stream/pitch,0);
%                     oscosa(k)=o_s(k)/cosd(job.alpha_rel_inlet);
%                     [rad_contr(k)] = mis_radial_contraction(directory);
%                     AtA1(k) = oscosa(k).*rad_contr(k);
                    
                    bl_inc = 1;
                    mix = mis_mix_out_average(mis_directory,bl_inc,0);
                    AtA1(k) = mix.AtA1;
                    AtA1_rhoV(k) = mix.AtA1_rhoV;
                    AtA1_P(k) = mix.AtA1_P;
                    % plot mises data summary
                    [p,h] = mis_plot_section(directory,h,color(i,:),plot_stuff,job);
                    % store Peak Mach number and loss (overall and freestream)
                    M_max(k) = p.M_max_ss;
                    Yp(k) = p.loss;
                    Yp_fr(k) = p.loss_fr;
                    dPo_Po(k)=Idat.dPo2_Po_mass;
                    
                end
                if i ~= length(perc_cam)
                k = k+1;
                continue
                end
            end
        end
        cd(directory)
        
        job.Axt = AxtAx1_sweep(m);

        % ADJUST FRACTION chord to betweeen covered passage and shock
        % (mainly for pitch-to-chord)
        [job]=ts_frac_chord_mises(job,[],2); % mode 2 - between shock and covered passage
        [q,h,job] = mis_run_section(directory,[],[0,0,1],0,1,0,job);
        if q.manual == 1
            disp(['****** Cannot run datum at : ' 'AxtAx1' num2str(round(100*AxtAx1_sweep(m))) 'PC' num2str(round(100*perc_cam(i))) ' ******'])
        end
        save([directory 'job.mat'],'job')
        
        Legend = [Legend ; {['AxtAx1 = ' num2str(round(100*AxtAx1_sweep(m)))]}];
        if plot_stuff == 1 && exist('AtA1','var')==1
            % PLOT FIGURES
            [h1,h2,h3]=plot_fig(h1,h2,h3,J,m,Legend);
            % clear all recalculated variables
            clear Yp Yp_fr M_max AtA1
            h = [];
        end

    end
end

%% sweep through chi_le

J = J + 1;
k = 1;
job = job_datum;
chi_le_datum = job_datum.chi_le;
dchi_le = [-1.5 -1 -0.5 0.5 1 1.5];
chi_le_sweep = chi_le_datum+dchi_le;
for m = 1:length(chi_le_sweep)
    
    % sweep through different perecentage cambers
    perc_cam_datum = job_datum.perc_cam;
    dperc_cam = [-0.2 -0.15 -0.1 -0.05 0 0.05 0.1 0.15 0.20];
    perc_min = 0; perc_max = 0.70;
    perc_cam = perc_cam_datum + dperc_cam;
    perc_cam = perc_cam(perc_cam>perc_min); perc_cam = perc_cam(perc_cam<perc_max);
    tot_cam_datum = job_datum.chi_le - job_datum.chi_te;
    
    for i = 1:length(perc_cam)

        job.chi_le = chi_le_sweep(m);

        job.perc_cam = perc_cam(i).*(job_datum.chi_le-job_datum.chi_te)./(job.chi_le-job.chi_te);
        directory = [job.new_datum_directory file_prefix '-' 'LE' num2str(round(100*chi_le_sweep(m))/100) 'PC' num2str(round(100*job.perc_cam)) '/MISES/'];
        job.rjm_directory = directory;
        % Make directory if required
        if exist(directory,'dir') == 0
            mkdir(directory);
            disp([directory ': does not exist'])
        else
            [Polarx, Ises] = mis_read_polarx('mises',directory);
            if isstruct(Polarx) ~= 1
                disp([directory ' : solution did not converge'])
                continue
            else
                if plot_stuff == 1
                    % calculate Athroat
                    plot_blade = 0; ote = 0; dev = []; cam_spline = 1;
                    [xrt,xrt_cam,xrt_ps,xrt_ss]=bl_construct_section(job,plot_blade,ote,dev,cam_spline);
                    mca.xrt_ss = xrt_ss;
                    mca.xrt_ps = xrt_ps;
                    pitch = 2*pi*job.r_mid/job.N; mt_stag_stream = [];
                    [~,o_s(k),xrt_throat]=ts_throat_calc(mca.xrt_ss,mca.xrt_ps,pitch,mt_stag_stream/pitch,0);
                    oscosa(k)=o_s(k)/cosd(job.alpha_rel_inlet);
                    [rad_contr(k)] = mis_radial_contraction(directory);
                    AtA1(k) = oscosa(k).*rad_contr(k);
                    % plot mises data summary
                    [p,h] = mis_plot_section(directory,h,color(i,:),plot_stuff,job);
                    title([file_prefix '-' 'TCPC'])
                    % store Peak Mach number and loss (overall and freestream)
                    M_max(k) = p.M_max_ss;
                    Yp(k) = p.loss;
                    Yp_fr(k) = p.loss_fr;
                end
                if i ~= length(perc_cam)
                    k = k+1;
                    continue
                end

            end
        end
        cd(directory)

        % ADJUST FRACTION chord to betweeen covered passage and shock
        % (mainly for pitch-to-chord)
        [job]=ts_frac_chord_mises(job,[],2); % mode 1 - between shock and covered passage
        [q,h,job] = mis_run_section(directory,[],[0,0,1],0,1,0,job);
        if q.manual == 1
            disp(['****** Cannot run datum at : ' 'LE' num2str(round(100*chi_le_sweep(m))/100) 'PC' num2str(round(100*job.perc_cam)) ' ******'])
        end
        save([directory 'job.mat'],'job')
        
        Legend = [Legend ; {['\chi_{le} = ' num2str(round(100*chi_le_sweep(m)))]}];
        if plot_stuff == 1 && exist('AtA1','var')==1
            % PLOT FIGURES
            [h1,h2,h3]=plot_fig(h1,h2,h3,J,m,Legend);
            % clear all recalculated variables
            clear Yp Yp_fr M_max AtA1
            h = [];
        end

    end
end

%% sweep through perc_cam

J = J + 1;

%% sweep through chi_te

J = J + 1;
k = 1;
job = job_datum;
chi_te_datum = job_datum.chi_te;
dchi_te = [-2 -1.5 -1 -0.5 0.5 1 1.5 2];
chi_te_sweep = chi_te_datum+dchi_te;
for m = 1:length(chi_te_sweep)
    
    % sweep through different perecentage cambers
    perc_cam_datum = job_datum.perc_cam;
    dperc_cam = [-0.2 -0.15 -0.1 -0.05 0 0.05 0.1 0.15 0.20];
    perc_min = 0; perc_max = 0.70;
    perc_cam = perc_cam_datum + dperc_cam;
    perc_cam = perc_cam(perc_cam>perc_min); perc_cam = perc_cam(perc_cam<perc_max);
    tot_cam_datum = job_datum.chi_le - job_datum.chi_te;
    
    for i = 1:length(perc_cam)

        job.chi_te = chi_te_sweep(m);

        job.perc_cam = perc_cam(i).*(job_datum.chi_le-job_datum.chi_te)./(job.chi_le-job.chi_te);
        directory = [job.new_datum_directory file_prefix '-' 'TE' num2str(round(100*chi_te_sweep(m)/100)) 'PC' num2str(round(100*job.perc_cam)) '/MISES/']
        job.rjm_directory = directory;
        % Make directory if required
        if exist(directory,'dir') == 0
            mkdir(directory);
            disp([directory ': does not exist'])
        else
            [Polarx, Ises] = mis_read_polarx('mises',directory);
            if isstruct(Polarx) ~= 1
                disp([directory ' : solution did not converge'])
                continue
            else
                if plot_stuff == 1
                    % calculate Athroat
                    plot_blade = 0; ote = 0; dev = []; cam_spline = 1;
                    [xrt,xrt_cam,xrt_ps,xrt_ss]=bl_construct_section(job,plot_blade,ote,dev,cam_spline);
                    mca.xrt_ss = xrt_ss;
                    mca.xrt_ps = xrt_ps;
                    pitch = 2*pi*job.r_mid/job.N; mt_stag_stream = [];
                    [~,o_s(k),xrt_throat]=ts_throat_calc(mca.xrt_ss,mca.xrt_ps,pitch,mt_stag_stream/pitch,0);
                    oscosa(k)=o_s(k)/cosd(job.alpha_rel_inlet);
                    [rad_contr(k)] = mis_radial_contraction(directory);
                    AtA1(k) = oscosa(k).*rad_contr(k);
                    % plot mises data summary
                    [p,h] = mis_plot_section(directory,h,color(i,:),plot_stuff,job);
                    title([file_prefix '-' 'TCPC'])
                    % store Peak Mach number and loss (overall and freestream)
                    M_max(k) = p.M_max_ss;
                    Yp(k) = p.loss;
                    Yp_fr(k) = p.loss_fr;
                end
                if i ~= length(perc_cam)
                    k = k+1;
                    continue
                end

            end
        end
        cd(directory)

        % ADJUST FRACTION chord to betweeen covered passage and shock
        % (mainly for pitch-to-chord)
        [job]=ts_frac_chord_mises(job,[],2); % mode 1 - between shock and covered passage
        [q,h,job] = mis_run_section(directory,[],[0,0,1],0,1,0,job);
        if q.manual == 1
            disp(['****** Cannot run datum at : ' 'TE' num2str(round(100*chi_te_sweep(m)/100)) 'PC' num2str(round(100*job.perc_cam)) ' ******'])
        end
        save([directory 'job.mat'],'job')
        
        Legend = [Legend ; {['\chi_{te} = ' num2str(round(100*chi_te_sweep(m)))]}];
        if plot_stuff == 1 && exist('AtA1','var')==1
            % PLOT FIGURES
            [h1,h2,h3]=plot_fig(h1,h2,h3,J,m,Legend);
            % clear all recalculated variables
            clear Yp Yp_fr M_max AtA1
            h = [];
        end

    end

end


    function [h1,h2,h3]=plot_fig(h1,h2,h3,J,m,Legend)
        figure(h1); grid on; hold on;
        plot(AtA1,M_max,'Marker',marker{J},'Markersize',8,'Linestyle','none','Color',color(m,:));
        %             plot(1-o_scosx,glob_inc,'Marker',marker{m},'Markersize',8,'Linestyle','none','Color',color(m,:),'MarkerFaceColor',color(m,:));
        ylabel('Peak Mach number')
        xlabel('A_{throat}/A_{inlet}')
        legend(Legend,'Location','Best','Fontsize',10)
        xlim([0.95,1.1])
        
        figure(h2); grid on; hold on;
        plot(AtA1,Yp,'Marker',marker{J},'Markersize',8,'Linestyle','none','Color',color(m,:));
        %             plot(1-o_scosx,glob_inc,'Marker',marker{m},'Markersize',8,'Linestyle','none','Color',color(m,:),'MarkerFaceColor',color(m,:));
        ylabel('Profile loss')
        xlabel('A_{throat}/A_{inlet}')
        legend(Legend,'Location','Best','Fontsize',10)
        xlim([0.95,1.1])
        
        figure(h3); grid on; hold on;
        plot(AtA1,Yp_fr,'Marker',marker{J},'Markersize',8,'Linestyle','none','Color',color(m,:));
        %             plot(1-o_scosx,glob_inc,'Marker',marker{m},'Markersize',8,'Linestyle','none','Color',color(m,:),'MarkerFaceColor',color(m,:));
        ylabel('Freestream loss')
        xlabel('A_{throat}/A_{inlet}')
        legend(Legend,'Location','Best','Fontsize',10)
        xlim([0.95,1.1])
    end
end