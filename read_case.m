function rcase = read_case(casename, type, run)

    if nargin < 2
        type = 'cpu';
    end

    if nargin < 3 || isempty(run)
        runpath = [];
    else
        runpath = ['run' num2str(run(end))];
    end

    base = pwd;

    switch type
    
        case 'cpu'
    
            fpath = fullfile(base, casename, runpath, 'input_cpu.txt');
            f = fopen(fpath);
        
            NB = str2num(fgetl(f));
            nprocs=0;
            npoints=0;
            for ib=1:NB
                nijk = str2num(char(split(fgetl(f))));
                nijk_procs = str2num(char(split(fgetl(f))));
                ptch = str2num(char(split(fgetl(f))));
                
                ni = nijk(1); nj = nijk(2); nk = nijk(3);
                
                npoints = npoints + ni*nj*nk;
        
                niproc = nijk_procs(1);
                njproc = nijk_procs(2);
                nkproc = nijk_procs(3);
        
                nprocs = nprocs + niproc*njproc*nkproc;
        
                im = ptch(1);
                ip = ptch(2);
                jm = ptch(3);
                jp = ptch(4);
        
                next_block{ib}.im = 0;
                next_block{ib}.ip = 0;
                next_block{ib}.jm = 0;
                next_block{ib}.jp = 0;
                
                next_patch{ib}.im = im;
                next_patch{ib}.ip = ip;
                next_patch{ib}.jm = jm;
                next_patch{ib}.jp = jp;
        
                if im == 0
                    temp = str2num(char(split(fgetl(f))));
                    next_block{ib}.im = temp(1);
                    next_patch{ib}.im = temp(2);
                end
                if ip == 0
                    temp = str2num(char(split(fgetl(f))));
                    next_block{ib}.ip = temp(1);
                    next_patch{ib}.ip = temp(2);
                end
                if jm == 0
                    temp = str2num(char(split(fgetl(f))));
                    next_block{ib}.jm = temp(1);
                    next_patch{ib}.jm = temp(2);
                end
                if jp == 0
                    temp = str2num(char(split(fgetl(f))));
                    next_block{ib}.jp = temp(1);
                    next_patch{ib}.jp = temp(2);
                end
        
                NI{ib} = ni;
                NJ{ib} = nj;
        
            end
        
            % Read corners
            temp = str2num(char(split(fgetl(f))));
            ncorner = temp(1);
            corner = {};
        
            for n=1:ncorner
                cor_type(n) = 0;
                temp = str2num(char(split(fgetl(f))));
                corner{n}.Nb = temp(1);
                cor_type(n) = temp(2);
        
                for m = 1:corner{n}.Nb
                    temp = str2num(char(split(fgetl(f))));
                    ib = temp(1);
                    ic = temp(2);
                    jc = temp(3);
        
                    corner{n}.block{m} = ib;
                    corner{n}.i{m} = ic;
                    corner{n}.j{m} = jc;
        
                end
        
            end
            
            solver.nk = nk;
        
            temp = str2num(char(split(fgetl(f))));
            solver.niter = temp(1);
            solver.nwrite = temp(2);
            solver.ncut = temp(3);
        
            temp = str2num(char(split(fgetl(f))));
            solver.cfl = temp(1);
            solver.sigma = temp(2);
            solver.ifsplit = temp(3);
            solver.ifsat = temp(4);
            solver.ifLES = temp(5);
        
            temp = str2num(char(split(fgetl(f))));
            bcs.Toin = temp(1);
            bcs.Poin = temp(2);
            bcs.pexit = temp(3);
            bcs.vin = temp(4);
            bcs.alpha = temp(5);
            bcs.cax = temp(6);
            bcs.aturb = temp(7);
            bcs.lturb = temp(8);
            bcs.ilength = temp(9);
            bcs.radprof = temp(10);
            bcs.gamma = 0.0;
            bcs.g_z = 0.0;
        
            % Gas props
            temp = str2num(char(split(fgetl(f))));
            gas.gam = temp(1);
            gas.cp = temp(2);
            gas.mu_ref = temp(3);
            gas.mu_tref = temp(4);
            gas.mu_cref = temp(5);
            gas.pr = temp(6);
        
            % Span, expan factor
            temp = str2num(char(split(fgetl(f))));
            solver.span = temp(1);
            blk.span = temp(1);
            solver.fexpan = temp(2);
        
            temp = str2num(char(split(fgetl(f))));
            solver.irestart = temp(1);
            solver.istats = temp(2);
        
            temp = str2num(char(split(fgetl(f))));
            n_in_blks = str2num(fgetl(f));
        
            for i=1:n_in_blks
                inlet_blks(i) = str2num(fgetl(f));
            end
        
            % Input for stability analysis
            temp = str2num(char(split(fgetl(f))));
            solver.istability = temp(1);
        
            % Input for incoming boundary layer and adiabatic/isothermal wall
            temp = str2num(char(split(fgetl(f))));
            bcs.twall = temp(3);
        
            fclose(f);
        
            if nk == 1
                solver.npp = ceil((npoints/nprocs)^(1/2));
            else
                solver.npp = ceil((npoints/nprocs)^(1/3));
            end
        
            blk = read_grid(casename);
            blk.span = solver.span;
            blk.npp = solver.npp;
        
            rcase.NB = NB;
            rcase.blk = blk;
            rcase.next_block = next_block;
            rcase.next_patch = next_patch;
            rcase.corner = corner;
            rcase.bcs = bcs;
            rcase.gas = gas;
            rcase.solver = solver;
            rcase.inlet_blocks = inlet_blks;

        case 'gpu'

            fpath = fullfile(base, casename, 'input_gpu.txt');
            f = fopen(fpath);
        
            tmp = str2num(char(split(fgetl(f))));
            NB = tmp(1); solver.nkprocs = tmp(2);
            nprocs=0;
            npoints=0;
            inlet_blks = [];

            for ib=1:NB
                nijk = str2num(char(split(fgetl(f))));
                ptch = str2num(char(split(fgetl(f))));
                
                ni = nijk(1); nj = nijk(2); nk = nijk(3);
                
                npoints = npoints + ni*nj*nk;
      
                im = ptch(1);
                ip = ptch(2);
                jm = ptch(3);
                jp = ptch(4);
        
                next_block{ib}.im = 0;
                next_block{ib}.ip = 0;
                next_block{ib}.jm = 0;
                next_block{ib}.jp = 0;
                
                next_patch{ib}.im = im;
                next_patch{ib}.ip = ip;
                next_patch{ib}.jm = jm;
                next_patch{ib}.jp = jp;
        
                if im == 0
                    temp = str2num(char(split(fgetl(f))));
                    next_block{ib}.im = temp(1);
                    next_patch{ib}.im = temp(2);
                end
                if ip == 0
                    temp = str2num(char(split(fgetl(f))));
                    next_block{ib}.ip = temp(1);
                    next_patch{ib}.ip = temp(2);
                end
                if jm == 0
                    temp = str2num(char(split(fgetl(f))));
                    next_block{ib}.jm = temp(1);
                    next_patch{ib}.jm = temp(2);
                end
                if jp == 0
                    temp = str2num(char(split(fgetl(f))));
                    next_block{ib}.jp = temp(1);
                    next_patch{ib}.jp = temp(2);
                end

                if im == 1
                    inlet_blks(end+1) = ib;
                end
        
                NI{ib} = ni;
                NJ{ib} = nj;
        
            end
        
            % Read corners
            temp = str2num(char(split(fgetl(f))));
            ncorner = temp(1);
            corner = {};
        
            for n=1:ncorner
                cor_type(n) = 0;
                temp = str2num(char(split(fgetl(f))));
                corner{n}.Nb = temp(1);
                %cor_type(n) = temp(2);
        
                for m = 1:corner{n}.Nb
                    temp = str2num(char(split(fgetl(f))));
                    ib = temp(1);
                    ic = temp(2);
                    jc = temp(3);
        
                    corner{n}.block{m} = ib;
                    corner{n}.i{m} = ic;
                    corner{n}.j{m} = jc;
        
                end
        
            end

            solver.nk = nk;

            % Read block groups

            nBlockGroups = str2num(char(fgetl(f)));


            for nbg = 1:nBlockGroups
                nb_bg(nbg) = str2num(char(fgetl(f)));
                temp = str2num(char(split(fgetl(f))));
                block{nbg} = temp;
            end

            
        
            temp = str2num(char(split(fgetl(f))));
            solver.niter = temp(1);
            solver.nwrite = temp(2);
            solver.ncut = temp(3);
        
            temp = str2num(char(split(fgetl(f))));
            solver.cfl = temp(1);
            solver.sigma = temp(2);
        
            temp = str2num(char(split(fgetl(f))));
            bcs.Toin = temp(1);
            bcs.Poin = temp(2);
            bcs.pexit = temp(3);
            bcs.vin = temp(4);
            bcs.alpha = temp(5);
            bcs.gamma = temp(6);
            bcs.aturb = temp(7);
            bcs.nturb = temp(8);
            bcs.iradprof = temp(9);
            bcs.g_z = temp(10);
        
            % Gas props
            temp = str2num(char(split(fgetl(f))));
            gas.gam = temp(1);
            gas.cp = temp(2);
            gas.mu_ref = temp(3);
            gas.mu_tref = temp(4);
            gas.mu_cref = temp(5);
            gas.pr = temp(6);
        
            % Span, expan factor
            temp = str2num(char(split(fgetl(f))));
            solver.span = temp(1);
            blk.span = temp(1);
            solver.fexpan = temp(2);
        
            temp = str2num(char(split(fgetl(f))));
            solver.irestart = temp(1);
            solver.istats = temp(2);
            solver.istability = 0;

            temp = str2num(char(split(fgetl(f))));
            solver.ilam = temp(1);
            bcs.theta = temp(2);
        
%             temp = str2num(char(split(fgetl(f))));
%             n_in_blks = str2num(fgetl(f));
%         
%             for i=1:n_in_blks
%                 inlet_blks(i) = str2num(fgetl(f));
%             end
%         
%             % Input for stability analysis
%             temp = str2num(char(split(fgetl(f))));
%             solver.istability = temp(1);
%         
%             % Input for incoming boundary layer and adiabatic/isothermal wall
%             temp = str2num(char(split(fgetl(f))));
%             bcs.twall = temp(3);
        
            fclose(f);
        
            if nk == 1
                solver.npp = ceil((npoints/nprocs)^(1/2));
            else
                solver.npp = ceil((npoints/nprocs)^(1/3));
            end
        
            blk = read_grid(casename);
            blk.span = solver.span;
            blk.npp = solver.npp;
            blk.nbg = nBlockGroups;
            blk.block_groups = block;
        
            rcase.NB = NB;
            rcase.blk = blk;
            rcase.next_block = next_block;
            rcase.next_patch = next_patch;
            rcase.corner = corner;
            rcase.bcs = bcs;
            rcase.gas = gas;
            rcase.solver = solver;
            rcase.inlet_blocks = inlet_blks;

    end

end



