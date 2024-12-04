% function to store AtA1, AtD, a1-x1 and Yp
% save as spline

clear all

%% get the correct filenames based on aero parameters

t_c_dat = 4/100;
s_c_dat = 0.80;


alpha_rel_inlet_dat = 60;
alpha_inlet_dat = 0;
M_inlet_dat = 0.95;
np = 0.95;
load_coef_dat = 0.40;

Ax_dat = 0.95;
% Axt_dat = 1.0;
dH_dat = 0.66;
mode = 2; % mode 1: constant psi & dh; mode 2: contant psi & Area ratio; mode 3: constnat dH and Area ratio ; mode 4: constant dH and throat area ratio
job.mode = mode;


sweep.t_c = [4 3 5]/100;
sweep.Axt = [0.975 1.025 1.0 0.95];
sweep.s_c = [0.60 0.70 0.80 0.90 1.0];
sweep.alpha_rel_inlet = [50 55 60 65];
sweep.Ax = [1.0 0.95 0.975 0.90];
sweep.load_coef = [0.40 0.35 0.30 0.45];
% extract variable that will be swept
var = fieldnames(sweep);

for k = 1:length(var)
    
    ext2 = [];
    filename = [];
    filestring = [];
    Legend = [];

    
    disp(['*** sweeping: ' var{k} ' ****'])
    
    for m = 1:length(sweep.(var{k}))
        
        
        clear Axt
        
        if strcmp(var{k},'s_c')==1
            s_c = sweep.s_c(m);
        else
            s_c = s_c_dat;
        end
        
        if strcmp(var{k},'t_c')==1
            t_c = sweep.t_c(m);
        else
            t_c = t_c_dat;
        end
        
        %% change inlet flow angle (or U/cpTo1) and Mach number
        if strcmp(var{k},'alpha_rel_inlet')==1
            alpha_rel = sweep.alpha_rel_inlet(m);
        else
            alpha_rel = alpha_rel_inlet_dat;
        end
        
        if strcmp(var{k},'alpha_inlet')==1
            alpha = sweep.alpha_inlet(m);
        else
            alpha = alpha_inlet_dat;
        end
        
        if strcmp(var{k},'M_inlet')==1
            Min = sweep.M_inlet(m);
        else
            Min = M_inlet_dat;
        end
        
        %% set polytropic efficiency
        job.np = np;
        
        %% change aeroparameters
        if strcmp(var{k},'load_coef')==1
            psi = sweep.load_coef(m);
        else
            psi = load_coef_dat;
        end
        if strcmp(var{k},'dH')==1
            dH = sweep.dH(m);
        else
            dH = dH_dat;
        end
        
        %% change streamtube contraction
        if strcmp(var{k},'Ax')==1
            Ax = sweep.Ax(m);
        else
            Ax = Ax_dat;
        end
        
        
        %% change streamtube contraction
        if strcmp(var{k},'Axt')==1
            Axt = sweep.Axt(m);
            %         else
            %             job.Axt = Axt;
        end
        
        if mode==1
            ext = ['_psi0' num2str(100*psi) 'A' num2str(round(100*Ax))];
        elseif mode  == 2
            if exist('Axt','var')==0
                ext = ['_psi0' num2str(100*psi) 'A' num2str(round(100*Ax))];
                ext_dat = ['_psi0' num2str(100*load_coef_dat) 'A' num2str(round(100*Ax_dat))];
            else
                ext = ['_psi0' num2str(100*psi) 'Axt' num2str(round(100*Axt)) 'A' num2str(round(100*Ax))];
                ext_dat = ['_psi0' num2str(100*load_coef_dat) 'A' num2str(round(100*Ax_dat))];
            end
        elseif mode == 3
            ext = ['_DH0' num2str(100*dH) 'A' num2str(round(100*Ax))];
        elseif mode == 4
            ext = ['_DH0' num2str(100*dH) 'Axt' num2str(round(100*Ax))];
        end
        %     storename = ['M' num2str(round(100*Min)) 'A' num2str(round(100*alpha_rel)/100) 'SC' num2str(round(100*s_c)) 'TC' num2str(round(100*t_c_sweep(m))) ext ext2]
        %     filestring = ['M' num2str(round(100*Min)) 'A' num2str(round(100*alpha_rel)/100) 'SC' num2str(round(100*s_c)) 'TC' num2str(round(100*t_c_sweep(m))) ext ext2]
        
        storename = ['M' num2str(round(100*Min)) 'A' num2str(round(100*alpha_rel)/100) 'SC' num2str(round(100*s_c)) 'TC' num2str(round(100*t_c)) ext ext2]
        filestring = ['M' num2str(round(100*Min)) 'A' num2str(round(100*alpha_rel)/100) 'SC' num2str(round(100*s_c)) 'TC' num2str(round(100*t_c))];
        datumname = ['M' num2str(round(100*M_inlet_dat)) 'A' num2str(round(100*alpha_rel_inlet_dat)/100) 'SC' num2str(round(100*s_c_dat)) 'TC' num2str(round(100*t_c_dat)) ext_dat];
        
        Legend = [Legend ; {[filestring ext(2:end)]}];
        % Search directory for filenames of different percentage cambers
        A = [dir('/mnt/Disk1/dl467/Documents/Test_Cases/MISES_solutions')];
        F = cell(length(A),1); for n = 1:length(A); F{n} = A(n).name; end;
        F(cellfun(@isempty,strfind(F,[filestring])) == 1) = [];
        F(cellfun(@isempty,strfind(F,[ext])) == 1) = [];
        if isempty(ext2)==0
            F(cellfun(@isempty,strfind(F,[ext2])) == 1) = [];
        else
            F(cellfun(@isempty,strfind(F,['-'])) == 0) = [];
        end
        filename = [filename ; {F} ]
        
        if exist(['/mnt/Disk2/dl467/Documents/MISES_summary/M' num2str(round(100*Min))])==0
            mkdir(['/mnt/Disk2/dl467/Documents/MISES_summary/M' num2str(round(100*Min))])
        end
        
        if exist(['/mnt/Disk2/dl467/Documents/MISES_summary/M' num2str(round(100*Min)) '/' datumname])==0
            disp(['**** Creating datum directory: ' '/mnt/Disk2/dl467/Documents/MISES_summary/M' num2str(round(100*Min)) '/' datumname ])
            mkdir(['/mnt/Disk2/dl467/Documents/MISES_summary/M' num2str(round(100*Min)) '/' datumname])
        end
        directory_store = ['/mnt/Disk2/dl467/Documents/MISES_summary/M' num2str(round(100*Min)) '/' datumname '/' var{k} '/'];
        if exist(directory_store)==0
            mkdir(directory_store)
        end
        
        if exist([directory_store storename '.mat'])~=0
            disp(['***** ' storename ' exists *****'])
            continue
        end
        
        %% store variables
        if exist('mis','var')==1
            fields = fieldnames(mis);
            for  i=1:length(fields)
                clear(fields{i})
                clear mis
            end
        end
        
        nn = length(filename{m});
        
        filename{m}
        
        i = 1;
        % re-order directory based on s/c
        for l=1:nn
            
            directory = ['/mnt/Disk1/dl467/Documents/Test_Cases/MISES_solutions/' filename{m}{l} '/MISES/' ]
            
            cd(directory)
            mis_directory = strrep(directory,'TURBOSTREAM','MISES');
            load([mis_directory 'job.mat']);
            
            %% calculate pre-shock Mach number and loss
            inc = [0];
            for q = 1:length(inc)
                
                % Read in flow file
                if exist([directory 'polarx.mises'],'file') ~= 0
                    [Polarx, Ises] = mis_read_polarx('mises',directory);
                else
                    %% break if solution not found
                    disp([directory  ' is empty: inc = ' num2str(inc(q))])
                    if job.initial_sol == 0
                        disp(['Initial solution did not converge'])
                    end
                    continue
                end
                
                % Check the point is converged
                [p,h] = mis_plot_section(directory,[],[0,0,1],0,job);
                if isfield(Polarx,'binl') == 0 || isempty(p) == 1
                    disp('Run Not Converged')
                    disp([directory  ' is empty: inc = ' num2str(inc(q))])
                    if job.initial_sol == 0
                        disp(['Initial solution did not converge'])
                    end
                    continue
                end
                
                
                % Find the design incidence from the polar run
                if exist('inc','var')==0 || isempty(inc) == 1
                    [~,j] = min(abs(Polarx.binl-Ises.binl));
                else
                    [~,j] = min(abs(Polarx.binl-(Ises.binl+inc)));
                    disp(['Actual incidence is ' num2str(Polarx.binl(j)-Ises.binl) ])
                end
                
                % store directory
                dr{i} = mis_directory;
                
                % obtain boundary layer and calculate Freeman and Cumpsty
                if exist([mis_directory 'idat.mises_01'],'file') ~= 0
                    Idat = mis_read_idat('mises_01',mis_directory);
                else
                    Idat = mis_read_idat('mises',mis_directory);
                end
                load([mis_directory 'section.mat']);
                
                % calculate pitch
                pitch = 2*pi*c.r_le/c.N;
                
                %% record AtA1 and AtDelta
                if isfield(job,'rad_le')==0
                    [xrt,sl_cx,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
                    xrt_ss = mca.xrt_ss; xrt_ps = mca.xrt_ps;
                    %%
                else
                    ote =0; dev =[]; cam_spline = 1;
                    [xrt,xrt_cam,xrt_ps,xrt_ss,cam_param,thick_param]=bl_construct_section(job,0,ote,dev,cam_spline);
                end
                % calculate value of o_scosa
                [~,o_s(i),xrt_throat]=ts_throat_calc(xrt_ss,xrt_ps,pitch,[],0);
                oscosa(i)=o_s(i)/cosd(job.f.Alpha+inc(q));
                [rad_contr(i)] = mis_radial_contraction(mis_directory);
                %             AtA1(i) = oscosa(i).*rad_contr(i);
                
                bl_inc = 1;
                mix = mis_mix_out_average(mis_directory,bl_inc,0);
                AtA1_bl(i) = mix.AtA1;
                AtA1_rhoV_bl(i) = mix.AtA1_rhoV;
                AtA1_P_bl(i) = mix.AtA1_P;
                P_throat_area_bl(i) = mix.P_area_throat;
                P_area_bl(i) = mix.P_area;
                P_bl(i) = mix.P;
                P_is_bl(i) = mix.P_is;
                P_is_rhoV_bl(i) = mix.P_is_rhoV;
                AtA1_area_bl(i) = mix.AtA1_area;
                bl_inc = 0;
                mix = mis_mix_out_average(mis_directory,bl_inc,0);
                AtA1(i) = mix.AtA1;
                AtA1_area(i) = mix.AtA1_area;
                AtA1_rhoV(i) = mix.AtA1_rhoV;
                AtA1_P(i) = mix.AtA1_P;
                P_throat_area(i) = mix.P_area_throat;
                P_area(i) = mix.P_area;
                P(i) = mix.P;
                P_is(i) = mix.P_is;
                P_is_rhoV(i) = mix.P_is_rhoV;
                
                if isfield(job,'chi_le')==1
                    inlet_cam(i) = job.chi_le;
                else
                    inlet_cam(i) = job.inlet_cam;
                end
                Athroat_Delta(i) = o_s(i).*rad_contr(i)./cosd(inlet_cam(i));
                %% global incidence
                if isfield(job,'f') == 1
                    inc(i) = job.f.Alpha + inc(q) - inlet_cam(i);
                end
                
                %% record peak Mach number and location
                % calculate maximum Mach number
                Cp = Polarx.cp{j}(:,1); s = Polarx.s{j}(:,1) / Polarx.s{j}(Polarx.iteb(1),1);
                Mis = Polarx.mn{j}(:,1);
                q = s > 0.05 & s < 0.99; Cp = Cp(q); s = s(q); Mis = Mis(q);
                [~, I] = min(Cp); I = max(I,3);
                
                
                % Fit a polynomial through points near peak suction to find the maximum more accurately
                s_temp = linspace(s(I-2),s(I+2),1000);
                Q = polyval(polyfit(s(I-2:I+2),Cp(I-2:I+2),3),s_temp);
                [~,i_max] = min(Q); s_Cp_max= s_temp(i_max);
                M = polyval(polyfit(s(I-2:I+2),Mis(I-2:I+2),3),s_temp);
                [~,i_max] = min(M); s_max(i) = s_temp(i_max);
                M_max(i) = M(i_max);
                
                
                %% record loss
                Yp(i) = Polarx.omega(j);
                Yp_fr(i) = Polarx.omega(j)-Polarx.omegv(j);
                dPo_Po(i)=Idat.dPo2_Po_mass;
                
                %% increment renumeration
                i = i + 1
                
            end
            
        end
        
        %% save
        if exist('AtA1','var') ~= 0
            mis.AtA1 = AtA1;
            mis.AtA1_area = AtA1_area;
            mis.AtA1_area_bl = AtA1_area_bl;
            mis.AtA1_rhoV = AtA1_rhoV;
            mis.AtA1_P = AtA1_P;
            mis.AtA1_bl = AtA1_bl;
            mis.AtA1_rhoV_bl = AtA1_rhoV_bl;
            mis.AtA1_P_bl = AtA1_P_bl;
            mis.Yp = Yp;
            mis.Yp_fr = Yp_fr;
            mis.dPo_Po = dPo_Po;
            mis.AtD = Athroat_Delta;
            mis.rad_contr = rad_contr;
            mis.inc = inc;
            mis.M_max = M_max;
            mis.s_max = s_max;
            mis.P_throat_area_avg = P_throat_area;
            mis.P = P;
            mis.P_area = P_area;
            mis.P_is = P_is;
            mis.P_is_rhoV = P_is_rhoV;
            mis.P_throat_area_avg_bl = P_throat_area_bl;
            mis.P_bl = P_bl;
            mis.P_area_bl = P_area_bl;
            mis.P_is_bl = P_is_bl;
            mis.P_is_rhoV_bl = P_is_rhoV_bl;
            mis.dr = dr;
            
            
            save([directory_store storename '.mat'],'mis')
        else
            disp(['***** ' storename ' does not exist or no converged solutions ****'])
            
        end
        
        
    end
    
end
        


