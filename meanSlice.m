classdef meanSlice < aveSlice
    % MEANSLICE Contains the 2D (spanwise averaged) mean flow


    properties
        time;
        nMean;
        meanTime;
        Pr;             % Turbulence production
        diss;
        pbar;
        Tbar;
        rev_gen_x;
        rev_gen_y;
        irrev_gen;
        diss_T;
        span;
        nk;
        rous;
        rovs;
        rows;
        k;
        advK;
        e_s = [];
        e_phi = [];
        e_irrev = [];
        e_rev = [];
        e_Pr = [];
        e_diss = [];
        roUddUdd;
        roVddVdd;
        roWddWdd;
        roUddVdd;
        roUddWdd;
        roVddWdd;
    end

    properties (Dependent = true, Hidden = true)
        eta_Kol             % Kolmogorov scale
        cellSize_Kol;       % Cell size/ Kolmogorov scale
        mut_opt;
        tau_Re;
        tau_Re_an;          % Anisotropic componant of Re stress tensor
        tau_an_mag;         % Magnitude of anisotropic componant of Re stress
        tau_Re_an_mag;      
        tau_Re_Boussinesq;  % Reynolds stress tensor calculated with 
        omega_opt;
        omega_opt_cleaned;
        Rij;
        tau_Re_mag;     % RijSij / sqrt(SijSij)
    end

    methods
        function obj = meanSlice(blk, gas, bcs, casedir, casetype, ishere, iPS)
            if nargin < 7
                iPS =false;
            end
            obj@aveSlice(blk, gas, bcs, iPS);
            disp('Constructing meanSlice')

            if nargin > 0

                if nargin > 3 && ~isempty(casedir)
                    if ishere
                        basedir = casedir;
                    else
                        basedir = fullfile(casedir, 'mean_flo');
                    end
                    fullfile(basedir, 'mean_time.txt')
    
                    if exist(fullfile(basedir,'nstats.txt'),'file')
                        fid = fopen(fullfile(basedir,'nstats.txt'));
                        nstats = str2double(fgetl(fid));
                        statstype = str2double(fgetl(fid));
                        fclose(fid);
                    else
                        nstats = 17;
                        statstype = 1;
                    end
                    nstats;
                    if isnan(statstype)
                        statstype = 2;
                    end
    
                    fid = fopen(fullfile(basedir, 'mean_time.txt'),'r');
                    while ~feof(fid) % Use lastest mean files
                        temp=fgetl(fid);
                    end
    %                 temp = fgetl(fid);
                    fclose(fid);
                    temp = str2num(temp);
                    obj.nMean = temp(1);
                    obj.meanTime = temp(3);
                    for nb = 1:obj.NB
    
                        ni = blk.blockdims(nb,1);
                        nj = blk.blockdims(nb,2);
        
                        ro = zeros(ni,nj);
                        ru = zeros(ni,nj);
                        rv = zeros(ni,nj);
                        rw = zeros(ni,nj);
                        Et = zeros(ni,nj);
    
                        ro2 = zeros(ni,nj);
                        rou2 = zeros(ni,nj);
                        rov2 = zeros(ni,nj);
                        row2 = zeros(ni,nj);
    
                        rouv = zeros(ni,nj);
                        rouw = zeros(ni,nj);
                        rovw = zeros(ni,nj);
    
                        p2 = zeros(ni,nj);
                        p = zeros(ni,nj);
                        T = zeros(ni,nj);
                        rous = zeros(ni,nj);
                        rovs = zeros(ni,nj);
                        rows = zeros(ni,nj);
                        diss = zeros(ni,nj);
                        qx_T = zeros(ni,nj);
                        qy_T = zeros(ni,nj);
                        qz_T = zeros(ni,nj);
                        irrev_gen = zeros(ni,nj);
                        diss_T = zeros(ni,nj);
    
                        sz = size(ro2);
                        icount(1:prod(sz)) = 0;
    
                        switch casetype
                            case 'cpu'
                                flopath = fullfile(basedir,  ['mean2_' num2str(nb) '_' num2str(obj.nMean)]);
                                flofile = fopen(flopath,'r');
                                fullfile(basedir, ['mnod2_' num2str(nb) '_' num2str(obj.nMean)]);
                                nodfile = fopen(fullfile(basedir, ['mnod2_' num2str(nb) '_' num2str(obj.nMean)]),'r');
                                A = fread(flofile,inf,'float64');
            %                     fprintf('%d %d %d\n',nb, length(A), nstats)
                                A = reshape(A,nstats,length(A)/nstats);
                                
                                B = fread(nodfile,inf,'uint32');
                                B = reshape(B,3,length(B)/3);
                        
                                fclose(flofile);
                                fclose(nodfile);
                                for n=1:size(A,2)
            
                                    i = B(1,n);
                                    j = B(2,n);
            
                                    if icount(sub2ind(sz,i,j)) == 0
            
                                        icount(sub2ind(sz,i,j)) = 1;
            
                                        ro(i,j) = A(1,n)/obj.meanTime;
                                        ru(i,j) = A(2,n)/obj.meanTime;
                                        rv(i,j) = A(3,n)/obj.meanTime;
                                        rw(i,j) = A(4,n)/obj.meanTime;
                                        Et(i,j) = A(5,n)/obj.meanTime;
                
                                        ro2(i,j) = A(6,n)/obj.meanTime;
                                        rou2(i,j) = A(7,n)/obj.meanTime;
                                        rov2(i,j) = A(8,n)/obj.meanTime;
                                        row2(i,j) = A(9,n)/obj.meanTime;
                
                                        rouv(i,j) = A(10,n)/obj.meanTime;
                                        rouw(i,j) = A(11,n)/obj.meanTime;
                                        rovw(i,j) = A(12,n)/obj.meanTime;
                
                                        p2(i,j) = A(13,n)/obj.meanTime;
            
                                        if statstype == 1
                                            p(i,j) = A(14,n)/obj.meanTime;
                                            T(i,j) = A(15,n)/obj.meanTime;
                                            rous(i,j) = A(16,n)/obj.meanTime;
                                        elseif statstype == 2
                                            rous(i,j) = A(14,n)/obj.meanTime;
                                            rovs(i,j) = A(15,n)/obj.meanTime;
                                            rows(i,j) = A(16,n)/obj.meanTime;
                                        end
            
                                        diss(i,j) = A(17,n)/obj.meanTime;
                                        if nstats > 17
                                            qx_T(i,j) = A(18,n)/obj.meanTime;
                                            qy_T(i,j) = A(19,n)/obj.meanTime;
                                            qz_T(i,j) = A(20,n)/obj.meanTime;
                                            irrev_gen(i,j) = A(21,n)/obj.meanTime; % Pr/(T^2*mu*cp) * Σqi*qi
                                        end
            
                                        if nstats > 21
                                            diss_T(i,j) = A(22,n)/obj.meanTime;
                                        end
                                    end
                                end
    
                            case 'gpu'
                                nstats_prim = 12;
                                if statstype == 4
                                    nstats_prim = 11;
                                end
                                nstats_budg = nstats-nstats_prim;
                                flopath = fullfile(basedir,  ['mean2_' num2str(nb) '_' num2str(obj.nMean)]);
                                flofile = fopen(flopath,'r');
                                fullfile(basedir, ['mnod2_' num2str(nb) '_' num2str(obj.nMean)]);
                                A = fread(flofile,ni*nj*nstats,'float64');
            %                     fprintf('%d %d %d\n',nb, length(A), nstats)
                                prim = reshape(A(1:ni*nj*nstats_prim),nstats_prim,[])';%length(A)/25)';
                                budg = reshape(A(ni*nj*nstats_prim+1:end),nstats_budg,[])';
                                
                                fclose(flofile);
    
                                nn = 1;
                                ro = reshape(prim(:,nn),ni,nj)/obj.meanTime;
                                nn = nn+1;
                                ru = reshape(prim(:,nn),ni,nj)/obj.meanTime;
                                nn = nn+1;
                                rv = reshape(prim(:,nn),ni,nj)/obj.meanTime;
                                nn = nn+1;
                                rw = reshape(prim(:,nn),ni,nj)/obj.meanTime;
                                nn = nn+1;
                                Et = reshape(prim(:,nn),ni,nj)/obj.meanTime;
        
                                if statstype ~= 4
                                    nn = nn+1;
                                    ro2 = reshape(prim(:,nn),ni,nj)/obj.meanTime;
                                end
                                nn = nn+1;
                                rou2 = reshape(prim(:,nn),ni,nj)/obj.meanTime;
                                nn = nn+1;
                                rov2 = reshape(prim(:,nn),ni,nj)/obj.meanTime;
                                nn = nn+1;
                                row2 = reshape(prim(:,nn),ni,nj)/obj.meanTime;
        
                                nn = nn+1;
                                rouv = reshape(prim(:,nn),ni,nj)/obj.meanTime;
                                nn = nn+1;
                                rouw = reshape(prim(:,nn),ni,nj)/obj.meanTime;
                                nn = nn+1;
                                rovw = reshape(prim(:,nn),ni,nj)/obj.meanTime;
        
                                nn = 1;
                                if statstype ~= 4
                                    p2 = reshape(budg(:,nn),ni,nj)/obj.meanTime;
                                    nn = nn+1;
                                end
    
                                rous = reshape(budg(:,nn),ni,nj)/obj.meanTime;
                                nn = nn+1;
                                rovs = reshape(budg(:,nn),ni,nj)/obj.meanTime;
                                if statstype ~= 4
                                    nn = nn+1;
                                    rows = reshape(budg(:,nn),ni,nj)/obj.meanTime;
                                end
                                
                                nn = nn+1;
                                diss = reshape(budg(:,nn),ni,nj)/obj.meanTime;
                                
                                nn = nn+1;
                                qx_T = reshape(budg(:,nn),ni,nj)/obj.meanTime;
                                nn = nn+1;
                                qy_T = reshape(budg(:,nn),ni,nj)/obj.meanTime;
                                nn = nn+1;
                                qz_T = reshape(budg(:,nn),ni,nj)/obj.meanTime;
                                
                                if statstype ~= 4
                                    nn = nn+1;
                                    irrev_gen = reshape(budg(:,nn),ni,nj)/obj.meanTime; % Pr/(T^2*mu*cp) * Σqi*qi
                                    nn = nn+1;
                                    diss_T = reshape(budg(:,nn),ni,nj)/obj.meanTime;
                                end
                            
                        end
    
                        Uf = ru./ro;
                        Vf = rv./ro;
                        Wf = rw./ro;
    
                        obj.ro{nb} = ro;
                        obj.u{nb} = Uf;
                        obj.v{nb} = Vf;
                        obj.w{nb} = Wf;
                        obj.Et{nb} = Et;
                        obj.pbar{nb} = (Et - 0.5*(rou2 + rov2 + row2))*(obj.gas.gam-1);
                        obj.Tbar{nb} = (obj.pbar{nb}.*obj.gas.gam)./(obj.gas.cp*(obj.gas.gam-1)*obj.ro{nb});
                        obj.diss{nb} = diss;
                        obj.rous{nb} = rous;
                        obj.rovs{nb} = rovs;
                        obj.rows{nb} = rows;
    
                        [DUDX,DUDY] = gradHO(blk.x{nb},blk.y{nb},obj.u{nb});
                        [DVDX,DVDY] = gradHO(blk.x{nb},blk.y{nb},obj.v{nb});
    
                        UddUdd = rou2./ro - Uf.*Uf;
                        VddVdd = rov2./ro - Vf.*Vf;
                        WddWdd = row2./ro - Wf.*Wf;
    
                        UddVdd = rouv./ro - Uf.*Vf;
                        UddWdd = rouw./ro - Uf.*Wf;
                        VddWdd = rovw./ro - Vf.*Wf;
    
                        obj.roUddUdd{nb} = rou2 - ro.*Uf.*Uf;
                        obj.roVddVdd{nb} = rov2 - ro.*Vf.*Vf;
                        obj.roWddWdd{nb} = row2 - ro.*Wf.*Wf;
    
                        obj.roUddVdd{nb} = rouv - ro.*Uf.*Vf;
                        obj.roUddWdd{nb} = rouw - ro.*Uf.*Wf;
                        obj.roVddWdd{nb} = rovw - ro.*Vf.*Wf;
    
                        obj.Pr{nb} = -ro.*(UddUdd.*DUDX + UddVdd.*(DUDY+DVDX) + VddVdd.*DVDY);
                        % obj.Pr{nb} = -obj.roUdd
                        obj.diss{nb} = diss;
                        obj.k{nb} = 0.5*(UddUdd + VddVdd + WddWdd); 
    
                        if nstats > 17
                            obj.rev_gen_x{nb} = qx_T;
                            obj.rev_gen_y{nb} = qy_T;
                            obj.irrev_gen{nb} = irrev_gen;
    
                        end
                        if nstats > 21
                            obj.diss_T{nb} = diss_T;
                        end
    
    
                    [drouk_dx,~]=gradHO(blk.x{nb},blk.y{nb},obj.ro{nb}.*obj.u{nb}.*obj.k{nb});
                    [~,drovk_dy]=gradHO(blk.x{nb},blk.y{nb},obj.ro{nb}.*obj.v{nb}.*obj.k{nb});
    
                    obj.advK{nb} = drouk_dx + drovk_dy;
                    end
                    
                    
                    obj.getBCs(blk.inlet_blocks{1});
                    
                end
                obj.span = blk.span;
                obj.nk = blk.nk;
            end
%             Mnow = obj.M;
%             Unow = obj.vel;
%             ronow = obj.ro;
%             munow = obj.mu;
%             %Mnow = Mnow{blk.inlet_blocks{1}};
%             pnow = obj.p;
%             %pnow = pnow{blk.inlet_blocks{1}};
%             
%             p0 = [];
%             Uinf = [];
%             muinf = [];
%             roinf = [];
%             for i=1:length(blk.inlet_blocks{1})
%                 p0now = pnow{blk.inlet_blocks{1}(i)}.*(1+((obj.gas.gam - 1)/2)*Mnow{blk.inlet_blocks{1}(i)}.^2).^(obj.gas.gam/(obj.gas.gam-1));
%                 p0 = [p0 p0now(40:100,:)];
%                 Uinf = [Uinf Unow{blk.inlet_blocks{1}(i)}(40:100,:)];
%                 muinf = [muinf munow{blk.inlet_blocks{1}(i)}(40:100,:)];
%                 roinf = [roinf ronow{blk.inlet_blocks{1}(i)}(40:100,:)];
%             end
%             obj.p0in = mean(p0,'all');
%             obj.Uinf = mean(Uinf,'all');
%             obj.muinf = mean(muinf,'all');
%             obj.roinf = mean(roinf,'all');
        end

        function addSlice(obj, newSlice)
            props2add = {"Pr", "diss", "rev_gen_x", "rev_gen_y", "irrev_gen", "diss_T", "p", "T", ...
                "rous", "rovs", "rows", "k", "advK", "ro", "u", "v", "w", "Et", ...
                "roUddUdd", "roVddVdd", "roWddWdd", "roUddVdd", "roUddWdd", "roVddWdd"};


            tot_time = obj.meanTime + newSlice.meanTime;

            for ip = 1:length(props2add)
                tmp = {};
                for ib = 1:obj.NB
                    if ~isempty(newSlice.(props2add{ip}))
                        tmp{ib} = (obj.(props2add{ip}){ib}*obj.meanTime + ...
                            newSlice.(props2add{ip}){ib}*newSlice.meanTime)/tot_time;
                    end
                end
                obj.(props2add{ip}) = tmp;
                
            end
            obj.meanTime = tot_time;
        end
        
        function value = get.eta_Kol(obj)
            nunow = obj.nu;
            dissnow = obj.diss;
            for ib = 1:obj.NB
                disstmp = max(1e-18*ones(size(dissnow{ib})), dissnow{ib});
                value{ib} = ((nunow{ib}.^3)./disstmp).^0.25;
            end
        end

        function value = get.cellSize_Kol(obj)
            etanow = obj.eta_Kol;
            csnow = obj.cellSize;
            for ib = 1:obj.NB
                value{ib} = csnow{ib}./etanow{ib};
            end
        end

        function value = get.mut_opt(obj)
            
            %Traceless strain tensor
            S = obj.St_an;
            % Traceless Re stress tensor
            tau = obj.tau_Re_an;
            for ib = 1:obj.NB

%                 [DUDX,DUDY] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.u{ib});
%                 [DVDX,DVDY] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.v{ib});

                
%                 S = zeros(obj.blk.blockdims(ib,1),obj.blk.blockdims(ib,2),3,3);
%                 tau = S;
%                 
%                 S(:,:,1,1) = 2*DUDX/3 - DVDY/3;
%                 S(:,:,2,2) = 2*DVDY/3 - DUDX/3;
%                 S(:,:,3,3) = -(DUDX+DVDY)/3;
% 
%                 S(:,:,1,2) = 0.5*(DUDY+DVDX);
%                 S(:,:,2,1) = S(:,:,1,2);

                % Traceless strain magnitude
%                 St = sqrt(sum(sum(S.*S,4),3));
                
%                 tau(:,:,1,1) = -2*obj.roUddUdd{ib}/3 + obj.roVddVdd{ib}/3 + obj.roWddWdd{ib}/3;
%                 tau(:,:,2,2) = -2*obj.roVddVdd{ib}/3 + obj.roUddWdd{ib}/3 + obj.roUddWdd{ib}/3;
%                 tau(:,:,3,3) = -2*obj.roWddWdd{ib}/3 + obj.roUddUdd{ib}/3 + obj.roVddVdd{ib}/3;
% 
%                 tau(:,:,1,2) = -obj.roUddVdd{ib};
%                 tau(:,:,2,1) = tau(:,:,1,2);

                num = sum(sum(tau{ib}.*S{ib},4),3);
                den = sum(sum(S{ib}.*S{ib},4),3);
                mask = obj.k{ib} < 10;
                num(mask) = 0;
                
%                 den  = max(den, 1000);
%                 num = obj.tau_Re{ib};
%                 den = obj.S_an_mag{ib};

                value{ib} = 0.5*abs(num./den);

            end
        end

        function value = mut_opt_ratio(obj)
            mut = obj.mut_opt;
            mu = obj.mu;
            for ib = 1:obj.NB
                value{ib} = mut{ib}./mu{ib};
            end
        end

        function value = get.tau_Re(obj)

            for ib = 1:obj.NB

                tau = zeros(obj.blk.blockdims(ib,1),obj.blk.blockdims(ib,2),3,3);
                
                tau(:,:,1,1) = -obj.roUddUdd{ib};
                tau(:,:,2,2) = -obj.roVddVdd{ib};
                tau(:,:,3,3) = -obj.roWddWdd{ib};

                tau(:,:,1,2) = -obj.roUddVdd{ib};
                tau(:,:,2,1) = tau(:,:,1,2);

                tau(:,:,1,3) = -obj.roUddWdd{ib};
                tau(:,:,3,1) = tau(:,:,1,3);

                tau(:,:,2,3) = -obj.roVddWdd{ib};
                tau(:,:,3,2) = tau(:,:,1,3);

                value{ib} = tau;
            end
        end

        function value = get.tau_Re_an(obj)

            tau_Re = obj.tau_Re;

            for ib = 1:obj.NB

                tr = tau_Re{ib}(:,:,1,1) + tau_Re{ib}(:,:,2,2) + tau_Re{ib}(:,:,3,3);

                tau = tau_Re{ib};
                for i = 1:3
                    tau(:,:,i,i) = tau(:,:,i,i) - tr/3;
                end

                value{ib} = tau;


%                 tau = zeros(obj.blk.blockdims(ib,1),obj.blk.blockdims(ib,2),3,3);
%                 
%                 tau(:,:,1,1) = -2*obj.roUddUdd{ib}/3 + obj.roVddVdd{ib}/3 + obj.roWddWdd{ib}/3;
%                 tau(:,:,2,2) = -2*obj.roVddVdd{ib}/3 + obj.roUddUdd{ib}/3 + obj.roWddWdd{ib}/3;
%                 tau(:,:,3,3) = -2*obj.roWddWdd{ib}/3 + obj.roUddUdd{ib}/3 + obj.roVddVdd{ib}/3;
% 
%                 tau(:,:,1,2) = -obj.roUddVdd{ib};
%                 tau(:,:,2,1) = tau(:,:,1,2);
% 
%                 value{ib} = sum(sum(abs(tau.*Snow{ib}),4),3);
            end
        end

        function value = get.tau_Re_an_mag(obj)

            tau = obj.tau_Re_an;
            for ib = 1:obj.NB
                value{ib} = sqrt(sum(sum(tau{ib}.*tau{ib},4),3));
            end

        end

        function value = get.Rij(obj)
            for ib = 1:obj.NB


                tau = zeros(obj.blk.blockdims(ib,1),obj.blk.blockdims(ib,2),3,3);
                
                tau(:,:,1,1) = obj.roUddUdd{ib};
                tau(:,:,2,2) = obj.roVddVdd{ib};
                tau(:,:,3,3) = obj.roWddWdd{ib};

                tau(:,:,1,2) = obj.roUddVdd{ib};
                tau(:,:,2,1) = tau(:,:,1,2);


                tau(:,:,1,3) = obj.roUddWdd{ib};
                tau(:,:,3,1) = tau(:,:,1,3);


                tau(:,:,2,3) = obj.roVddWdd{ib};
                tau(:,:,3,2) = tau(:,:,2,3);

                value{ib} = tau ./ obj.ro{ib};
            end
        end

        function value = get.tau_Re_mag(obj)
            R = obj.Rij;
            S = obj.St;
            for ib = 1:obj.NB

                value{ib} = sum(sum(R{ib}.*S{ib},4),3) ./ ...
                    sqrt(sum(sum(S{ib}.*S{ib},4),3));

            end
        end

        function value = get_ctau(obj)
            disp('Calculating c_tau from Rij and Sij')
            inds = obj.BLedgeInd;
            Uenow = obj.Ue;
            ronow = obj.oGridProp('ro');
            for i=1:length(inds)
                roe(i) = ronow(i,inds(i));
            end
            
            tau = obj.oGridProp('tau_Re_mag');
            value = -tau./(roe'.*Uenow.^2);
        end

        function write_mut_opt_input(obj, folder)

            mto = obj.mut_opt;
            for ib = 1:obj.NB

                tau_Re_tr = 2*obj.ro{ib}.*obj.k{ib}/3;
                
                n = 1;
                A = [];
                B = [];

                for i=1:obj.blk.blockdims(ib,1)
                    for j=1:obj.blk.blockdims(ib,2)
                        A(n) = mto{ib}(i,j);
                        A(n+1) = tau_Re_tr(i,j);
                        B(n) = i;
                        B(n+1) = j;
                        n = n+2;
                    end
                end
                
                mto_file = fopen(fullfile(folder, ['mut_opt_' num2str(ib)]),'w');
                nod_file = fopen(fullfile(folder, ['mto_nod_' num2str(ib)]),'w');
                
                fwrite(mto_file,A,'float64');
                fwrite(nod_file,B,'int32');
                fclose(mto_file);
                fclose(nod_file);
            end
        end


        function value = get.omega_opt(obj)
            mto = obj.mut_opt;
            mul = obj.mu;
            ronow = obj.ro;
            know = obj.k;

            for ib = 1:obj.NB
                om = ronow{ib}.*know{ib}./max(mul{ib},mto{ib});
                value{ib} = abs(om);
                
            end
        end

        function value = get.omega_opt_cleaned(obj)
            mto = obj.mut_opt_cleaned;
            mul = obj.mu;
            ronow = obj.ro;
            know = obj.k;

            del_all = obj.delta99;
            offset=1;

            for ib = 1:obj.NB
                
                del = del_all(offset:offset+obj.blk.blockdims(ib,1)-1);
                offset = offset+obj.blk.blockdims(ib,1)-1;

                om = ronow{ib}.*know{ib}./max(mul{ib},mto{ib});
                value{ib} = abs(om);

                ynow = obj.blk.y{ib};

                for i=1:obj.blk.blockdims(ib,1)
                    t = 1 - tanh(((ynow(i,:)/del(i) - 0.9)*10));
                    value{ib}(i,:) = t.*value{ib}(i,:);
                end

                                
            end

        end


        function value = get_p(obj)
            value = obj.pbar;
        end

        function value = get_T(obj)
            value = obj.Tbar;
        end

        function obj.set_p(obj,value)
            obj.pbar = value;
        end

        function obj.set_T(obj,value)
            obj.Tbar = value;
        end

        function value = get_mut(obj)
            value = obj.mut_opt;
        end
        
        function value = mut_opt_cleaned(obj, smoothwindow)

            if nargin < 2
                smoothwindow = 1;
            end
            del_all = obj.delta99;
            offset=1;

            for ib=1:obj.NB
                mtr = obj.mut_ratio{ib}; 
                ynow = obj.blk.y{ib};

                del = del_all(offset:offset+obj.blk.blockdims(ib,1)-1);
                offset = offset+obj.blk.blockdims(ib,1)-1;
    
                for i=1:obj.blk.blockdims(ib,1)
                    
                    mtr_max = max(mtr(i,ynow(i,:) < 0.9*del(i)));
                    lim = mtr_max * (1-tanh(((ynow(i,:)/del(i) - 0.9)*10)));
                    mtr(i,:) = min(mtr(i,:),lim);
                end
    
                if smoothwindow > 1
                    mtr = smoothdata2(mtr, 'lowess',smoothwindow);
                end
    
                value{ib} = mtr.*obj.mu{ib};
            end

        end

        function value = mtr_cleaned(obj, smoothwindow)

            if nargin < 2
                smoothwindow = 1;
            end
            del_all = obj.delta99;
            offset=1;

            for ib=1:obj.NB
                mtr = obj.mut_ratio{ib}; 
                ynow = obj.blk.y{ib};

                del = del_all(offset:offset+obj.blk.blockdims(ib,1)-1);
                offset = offset+obj.blk.blockdims(ib,1)-1;
    
                for i=1:obj.blk.blockdims(ib,1)
                    
                    mtr_max = max(mtr(i,ynow(i,:) < 0.9*del(i)));
                    lim = mtr_max * (1-tanh(((ynow(i,:)/del(i) - 0.9)*10)));
                    mtr(i,:) = min(mtr(i,:),lim);
                end
    
                if smoothwindow > 1
                    mtr = smoothdata2(mtr, 'lowess',smoothwindow);
                end
    
                value{ib} = mtr;
            end

        end

        
        function write_as_inst_slice(obj, path)

            for ib = 1:obj.NB

                f = fopen(fullfile(path, ['flow_' num2str(ib)]),'w');

                A(:,1) = reshape(obj.ro{ib},[],1);
                A(:,2) = reshape(obj.u{ib},[],1).*A(:,1);
                A(:,3) = reshape(obj.v{ib},[],1).*A(:,1);
                A(:,4) = reshape(obj.w{ib},[],1).*A(:,1);
                A(:,5) = reshape(obj.Et{ib},[],1);

                A = reshape(A',[],1);
                fwrite(f,A,'float64');
                fclose(f);

            end

        end

        function probes = xyPlusToProbeInds(obj, x, yp, k)

            if nargin < 4
                k=1;
            end

            probes = {};

            for i = 1:length(yp)
                [probe.i, probe.j, probe.nb] = grid_inds_at_y_plus(obj,x,yp(i));
                probe.k = k;
                probes{end+1} = probe;
            end

        end
        
    end


end
