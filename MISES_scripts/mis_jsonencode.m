% function to store varibales for Aljaz after mis_plot_summary & mis_store_AtA1 has been run

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
%% extract variable that will be swept
% var = fieldnames(sweep);
var = {'t_c','Axt','Ax'};
%%
nk = 1;
for k = 1:length(var)
    
    var{k}
    
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
            Axt = sweep.Axt(m)
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
        
        
        %% collect metadata
        directory_store = ['/mnt/Disk2/dl467/Documents/MISES_summary/M' num2str(round(100*Min)) '/' datumname '/' var{k} '/'];
        if exist(directory_store)==0
            disp(['****** Metadata directory does not exist yet for ' num2str(storename) ' ****'])
            continue
        end
        
        if exist([directory_store storename '.mat'])~=0
            disp(['***** ' storename ' exists *****'])
            load([directory_store storename '.mat']);
            fields = fieldnames(mis);
        else
            disp(['****** Metadata mat file does not exist yet for ' num2str(storename) ' ****'])
            continue
        end
       
        %% collect plot data

        % re-order directory based on s/c
        for l=1:length(mis.dr)
            
            directory = mis.dr{l};
            
            % rename mis.dr directory to just name
            name = strrep(mis.dr(l),'/mnt/Disk1/dl467/Documents/Test_Cases/MISES_solutions/',''); 
            name = strrep(name,'/MISES/',''); 
            mis_dr.name = name; 
            
            cd(directory)
            mis_directory = strrep(directory,'TURBOSTREAM','MISES');
            load([mis_directory 'job.mat']);
            
            for j = 1:length(fields)
                if l<=length(mis.(fields{j}))   
                mis_dr.(fields{j}) = mis.(fields{j})(l);
                end
            end

            % Read in flow file
            if exist([directory 'polarx.mises'],'file') ~= 0
                [Polarx, Ises] = mis_read_polarx('mises',directory);
            else
                disp('File Not Found')
                p = [];
            end
            
            % Check the point is converged
            if isfield(Polarx,'binl') == 0
                disp('Run Not Converged')
                p = [];
            end
            
            % Read in the grid coodinates
            if exist([directory 'idat.mises_01'],'file') ~= 0
                Idat = mis_read_idat('mises_01',directory);
            else
                Idat = mis_read_idat('mises',directory);
            end
            
            
            % Find the design incidence from the polar run
            [~,j] = min(abs(Polarx.binl-Ises.binl));
           
           
            %% get isentropiuc Mach number solutions
            
            % Record the boundary layer parameters
            q_1 = Polarx.ileb(1):Polarx.iteb(1);
            q_2 = Polarx.ileb(2):Polarx.iteb(2);
            s_1 = Polarx.s{j}(q_1,1);
            s_1 = s_1 / max(s_1);
            s_2 = Polarx.s{j}(q_2,2); s_2 = s_2 / max(s_2);
            M.s_1 = s_1; M.s_2 = s_2;
            
            Cp_1 = Polarx.cp{j}(q_1,1); Cp_2 = Polarx.cp{j}(q_2,2);
            Mis_1 = Polarx.mn{j}(q_1,1); Mis_2 = Polarx.mn{j}(q_2,2);
            M.Mis_1 = Mis_1; M.Mis_2 = Mis_2;
            
            %% get blade geometry
            load([job.rjm_directory 'section.mat']) % contains c.xrt and c.N and c.r_le and c.r_te
            r_mid = 0.5*(c.r_le + c.r_te);
            pitch = 2*pi*r_mid/c.N;
            p.xrt = c.xrt;
            p.xrt_pos_pitch(:,1) = c.xrt(:,1);
            p.xrt_neg_pitch(:,1) = c.xrt(:,1);
            p.xrt_pos_pitch(:,2) = c.xrt(:,2) + pitch;
            p.xrt_neg_pitch(:,2) = c.xrt(:,2) - pitch;
            
            %% get throat plane
            bl_inc = 1;
            %% if you want to calculate with boundary layer %%
            if exist('bl_inc','var')==0 || isempty(bl_inc) == 1 || bl_inc == 0
                bl=[];
            elseif bl_inc == 1
                bl = [Idat.x(:,1).*c.r_le Idat.y(:,1).*c.r_le+pitch];
                x = bl(:,1);
                y = bl(:,2);
                % construct section
                ote = 1; dev =[];
                cam_spline = 1;
                % construct section
                %% this is legacy
                if isfield(c,'xs_sb')==0
                    [xrt,sl_cx,mca]=mca_camber(job.inlet_cam,job.tot_cam,job.percentage_cam,job.pitch_chord,job.ss_points,job.chord,job.t,job.R_le,job.degree,job.knot_coef,0);
                    xrt_ss = mca.xrt_ss; xrt_ps = mca.xrt_ps;
                    %%
                else
                    [xrt,xrt_cam,xrt_ps,xrt_ss]=bl_construct_section(job,0,ote,dev,cam_spline);
                end
                
                y = y(x>min(xrt_ss(:,1)/pitch));
                x = x(x>min(xrt_ss(:,1)/pitch));
                y = y(x<max(1.05*xrt_ss(:,1)/pitch));
                x = x(x<max(1.05*xrt_ss(:,1)/pitch));
                bl = [x,y];
            end
            mode_throat = directory; plot_throat = 1; pitch_ref = []; xrt_stag_line = [Idat.x(:,end) Idat.y(:,end)].*c.r_le; alpha_in = abs(job.alpha_rel_inlet);
            ote =0; dev =[]; cam_spline = 1;
            [xrt,xrt_cam,xrt_ps,xrt_ss]=bl_construct_section(job,0,ote,dev,cam_spline);
            [~,throat_min,xrt_throat,xrt_planes] = ts_throat_calc(xrt_ss,xrt_ps,pitch,bl,plot_throat,pitch_ref,xrt_stag_line,alpha_in,mode_throat);
            p.xrt_throat_bl = xrt_throat.*pitch;
            p.xrt_stag_line = xrt_stag_line;
            p.bl = bl;
            
            %% get contours and blade geometry
   
            H=figure; grid on; hold on;
           
            H.Position = [471.0000 143 331.0000 603.7772];

            
            % determine contour levels
            
            M_contours = [0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.05 1.1 1.15 1.2 1.25 1.30 1.35 1.40];
            % plot rest of contours
            [C,f]=contour(Idat.x_cell.*c.r_le,Idat.y_cell.*c.r_le,Idat.M,M_contours,'Color',[0.25,0.25,0.25]);
            p.C = C;            
            [C_pitch,f]=contour(Idat.x_cell.*c.r_le,Idat.y_cell.*c.r_le+pitch,Idat.M,M_contours,'Color',[0.25,0.25,0.25]);
            p.C_pitch = C_pitch;
            close all
            
            %% get camber parameters and suction surface
            cam_spline = 1;
            cam_param = bl_construct_camber(job,[],cam_spline);
            cam.camber = 1-cam_param.cam;
            cam.s_cam = cam_param.s;
            
            mt = c.mt(c.mt(:,1) < c.mt_le_cen(1),:);
            [~,i_ss]=min(c.mt(:,1)-c.mt_le(1));
            % c.mt_ss = flip(c.mt(1:i_ss,:));
            c.mt_ss = c.mt(1:i_ss,:);
            c.mt_ps = c.mt(i_ss+1:end,:);
           
            mt_ss = c.mt(c.mt_ss(:,1) > c.mt_le_cen(1),:);
            mt_ss = flip(mt_ss);
            mt_ps = c.mt(c.mt_ps(:,1) > c.mt_le_cen(1),:);
            mt_ps = flip(mt_ps);
            
            scl_ss = [0 ; cumsum(sum(diff(mt_ss,1,1).^2,2).^0.5,1)];
            scl_ss = scl_ss/max(scl_ss);
            scl_ps = [0 ; cumsum(sum(diff(mt_ps,1,1).^2,2).^0.5,1)];
            scl_ps = scl_ps/max(scl_ps);
            
            
            % Calculate the angle the stagnation point makes with the leading edge
            dtdx = grad_mg(mt(:,1),mt(:,2));
            psi = atand(dtdx);
            dtdx_ss = grad_mg(mt_ss(:,1),mt_ss(:,2));
            psi_ss = atand(dtdx_ss);
            dtdx_ps = grad_mg(mt_ps(:,1),mt_ps(:,2));
            psi_ps = atand(dtdx_ps);
            cam.theta_ss = psi_ss;
            cam.s_ss = scl_ss;
            cam.theta_ps = [];
            cam.s_ps = [];
            
             %% get streamline
             % read stream file
             [stream] = mis_read_stream('mises',mis_directory);
            
            

            %% store vriables in class
            A{nk}.metadata = mis_dr;
            A{nk}.contour = p;
            A{nk}.distribution = M;
            A{nk}.camber = cam;
            A{nk}.stream = stream;
%             A{nk} = {mis_dr ; p ; M ;  cam_param };
            
            
            nk = nk +1;
            
        end
        
    end
    
end


%% save structure in jsoncode for Aljaz
Aljaz = jsonencode(A);
save(['/mnt/Disk2/dl467/Documents/dbslice/' datumname '_t_c_Axt_Ax.mat'],'Aljaz')
fid = fopen(['/mnt/Disk2/dl467/Documents/dbslice/' datumname '_t_c_Axt_Ax.json'],'w');
fprintf( fid, '%s', Aljaz )
fclose(fid)
cd(['/mnt/Disk2/dl467/Documents/dbslice/'])

        


