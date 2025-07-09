classdef DNS_case < handle
    %  DNS_CASE Class containing 3DNS case information and methods to read
    %  and postprocess results

    properties
        casename;
        casepath;
        casetype;
        topology;
        runpath;
        runpaths;
        run;
        blk;
        solver;
        bcs;
        gas;
        % x;
        % y;
        nSlices;
        nProbes;
        nSkip;
        probes;
        NB;
        y_inlet;
        inlet_width;
        nj_inlet;
        meanFlow = meanSlice.empty;
        meanFlowPS = meanSlice.empty;
        ransFlow = RANSSlice.empty;
        meanFlows = {};
        instFlow = volFlow.empty;
        startFlow = spanAveSlice.empty;
        endFlow = spanAveSlice.empty;
        iSlices;
        jSlices = jSlice.empty;
        kSlices = kSlice.empty;
        inflowTurb = volTurbulence.empty;
        spanAveFlow = kSlice.empty
        inletProf = inletProfile.empty;
        RANSSlices;
        trip;
        iTrip = false;
        e_unst = [];
        cell_area = [];
        if_rans;
        speed_up;
        pitch;
    end

    properties (Dependent = true, Hidden=true)
        ftrip;
        Re_k;        % Trip Reynolds number
        trip_Re_x;
        aspect_ratio;
    end

    methods
        function obj = DNS_case(casename,run,args)
            %DNS_CASE Construct an instance of this class
            %   Detailed explanation goes here
            if nargin > 0 && ~isempty(casename)
                tmp = split(casename, '/');
                if length(tmp) > 1
                    obj.casename = tmp{end};
                    obj.casepath = casename;
                else
                    obj.casename = casename;
                    obj.casepath = fullfile(pwd,obj.casename);
                end
                obj.runpaths = {};
                
                if nargin < 2 || isempty(run)
                    obj.runpath = obj.casepath;
                    obj.runpaths{1} = obj.runpath;
                    obj.run = [];
                elseif length(run) == 1
                    obj.run = run;
                    obj.runpath = fullfile(obj.casepath,['run' num2str(obj.run)]);
                    obj.runpaths{1} = obj.runpath;
                else
                    obj.run = run;
                    obj.runpath = fullfile(obj.casepath,['run' num2str(obj.run(end))]); 
                    for ir = 1:length(run)
                        obj.runpaths{ir} = fullfile(obj.casepath,['run' num2str(obj.run(ir))]);
                    end
                end
                disp(obj.casepath)

                if exist(fullfile(obj.casepath,'body.txt'),'file')
                    obj.casetype = 'gpu';
                else
                    obj.casetype = 'cpu';
                end
                fprintf('Case type: %s\n', obj.casetype);


                if exist(fullfile(obj.casepath,'rans.txt'),'file')
                    f = fopen(fullfile(obj.casepath,'rans.txt'),'r');
                    try
                        obj.if_rans = str2num(fgetl(f));
                        obj.speed_up = str2num(fgetl(f));
                    catch
                        obj.if_rans = 0;
                        obj.speed_up = 0;
                    end
                    fclose(f)
                else
                    obj.if_rans = false;
                    obj.speed_up = 1.0;
                end

                rcase = read_case(casename, obj.casetype, obj.run);
                obj.NB = rcase.NB;
                obj.blk = rcase.blk;
                obj.blk.next_block = rcase.next_block;
                obj.blk.next_patch = rcase.next_patch;
%                 obj.next_block = rcase.next_block;
%                 obj.next_patch = rcase.next_patch;
                obj.blk.corner = rcase.corner;
                obj.bcs = rcase.bcs;
                obj.gas = rcase.gas;
                obj.gas.rgas = obj.gas.cp * (obj.gas.gam - 1)/obj.gas.gam;
                obj.solver = rcase.solver;
                obj.blk.inlet_blocks{1} = rcase.inlet_blocks;
                obj.blk.outlet_blocks{1} = rcase.outlet_blocks;
                obj.blk.z = linspace(0, obj.solver.span, obj.blk.nk);
                obj.blk.viewarea = [];
                obj.pitch = 0;
                if nargin < 3 || ~isfield(args,'topology')
                    if obj.NB == 9
                        obj.topology = 1;
                    elseif obj.NB == 12
                        obj.topology = 2;
                        obj.blk.viewarea = [-0.6 2 -0.5 0.5];
                    end
                else
                    obj.topology = args.topology;
                end

                obj.construct_o_blocks;
                
                if ~isempty(obj.blk.viewarea)
                    obj.blk.aspect = [(obj.blk.viewarea(2)-obj.blk.viewarea(1)) ...
                        (obj.blk.viewarea(4)-obj.blk.viewarea(3)) 1];
                end
    
                if isfile(fullfile(obj.runpath,'slice_time.txt'))
                    obj.nSlices = size(readmatrix(fullfile(obj.runpath,'slice_time.txt')),1);
                end
    
                obj.solver.nk = obj.blk.nk;
    
                if length(obj.blk.inlet_blocks) == 1
                    obj.nj_inlet = 0;
                    y_inlet = [];
                    for i=1:length(obj.blk.inlet_blocks{1})
                        y_inlet = [y_inlet obj.blk.y{obj.blk.inlet_blocks{1}(i)}(1,:,1)];
                        ymidnow(i) = obj.blk.y{obj.blk.inlet_blocks{1}(i)}(1,ceil(end/2));
                    end
                    obj.y_inlet = sort(unique(y_inlet));
                    obj.nj_inlet = length(obj.y_inlet);
                    obj.inlet_width = max(obj.y_inlet) - min(obj.y_inlet);
                    
                    [~,inds] = sort(ymidnow);
                    obj.blk.inlet_blocks{1} = obj.blk.inlet_blocks{1}(inds);
                end
                
    
                obj.trip = [];
                if obj.solver.istability == 2
                    obj.iTrip = true;
                    obj.readTrip;
                else
                    obj.iTrip = false;
                end
                slices = dir(fullfile(obj.runpath,'k_cuts','kcu2_1_*'));
                obj.nSlices = length(slices);
                if obj.NB == 12
                    obj.blk.io_surfaces.blks = [1 2 2 7 11 11 10 12 12 8 3 3];
                    obj.blk.io_surfaces.types = [1 1 3 3 3 2 2 2 4 4 4 1];
                end
            end
        end

        function writeCase(obj)
            obj.writeGridFiles;
            obj.writeInputFiles;
            if ~isempty(obj.instFlow)
                obj.instFlow.writeFlow();
            end
        end

        function construct_o_blocks(obj)
            % Find first wall block, assume flow im -> ip
            for ib = 1:obj.NB
                if (obj.blk.next_block{ib}.jm == 0 ...
                        && obj.blk.next_patch{ib}.jm == 3) || ...
                    (obj.blk.next_block{ib}.jp == 0 ...
                        && obj.blk.next_patch{ib}.jp == 3)
                    obj.blk.oblocks = ib;
                    obj.blk.oblocks_flip = 0;
                    break
                end
            end

            % Construct wall blocks list and connectivity
            nb = obj.blk.next_block{ib}.ip;
            np = 'ip';
            while nb ~= 0 && ~ismember(nb, obj.blk.oblocks)
                obj.blk.oblocks(end+1) = nb;
                if obj.blk.next_patch{ib}.(np) == 1
                    obj.blk.oblocks_flip(end+1) = 0;
                    nb = obj.blk.next_block{nb}.ip;
                    np = 'ip';
                else
                    obj.blk.oblocks_flip(end+1) = 1;
                    nb = obj.blk.next_block{nb}.im;
                    np = 'im';
                end
                ib = obj.blk.oblocks(end);
            end
            
        end

        function calculate_wall_distance(obj)
            ib = obj.blk.oblocks(1);
            if (obj.blk.next_block{ib}.jp == 0 ...
                        && obj.blk.next_patch{ib}.jp == 3)
                jsurf = size(obj.blk.x{ib},2);
            else
                jsurf = 1;
            end
            xprof = [];
            yprof = [];
            for ib = obj.blk.oblocks
                xprof = [xprof obj.blk.x{ib}(:,jsurf)'];
                yprof = [yprof obj.blk.y{ib}(:,jsurf)'];
            end
            if obj.pitch ~= 0
                xprof = [xprof xprof];
                yprof = [(yprof - obj.pitch) yprof];
            end
            xnow = obj.blk.x;
            ynow = obj.blk.y;
            [~,mb] = max(prod(obj.blk.blockdims(:,1:2),2));
            parfor ib = 1:obj.NB
                [ni, nj] = size(xnow{ib});
                for i=1:ni
                    if ib == mb && mod(i,20) == 0
                        fprintf('i=%d/%d\n',[i,ni])
                    end
                    for j=1:nj
                        dist = sqrt((xprof - xnow{ib}(i,j)).^2 + (yprof - ynow{ib}(i,j)).^2);
                        walldist{ib}(i,j) = min(dist);
                    end
                end
            end
            obj.blk.walldist = walldist;
        end

        function slice = readSingleKSlice(obj,numslice)
            switch obj.casetype
                case 'cpu'
                    slicetime = readmatrix(fullfile(obj.runpath,'slice_time.txt'));
                case 'gpu'
                    slicetime = readmatrix(fullfile(obj.runpath,'kslice_time.txt'));
            end
            slicenums = slicetime(end-obj.nSlices+1:end,1);
            slicenum = slicenums(numslice);
            slice =kSlice(obj.blk,obj.gas,obj.bcs,obj.runpath,slicenum);
        end

        function slice = readSingleJSlice(obj,numslice)
            slicetime = readmatrix(fullfile(obj.runpath,'slice_time.txt'));
            slicenums = slicetime(end-obj.nSlices+1:end,1);
            slicenum = slicenums(numslice);
            slice =jSlice(obj.runpath,slicenum,obj.blk,obj.gas);
        end

        function readSpanAveFlo(obj)
            obj.spanAveFlow = kSlice(obj.blk,obj.gas,obj.bcs,obj.runpath,'flow_2d',[],obj.casetype,true);
        end

        function readRansFlow(obj)
            obj.ransFlow = RANSSlice(obj.blk, obj.gas, obj.bcs, '3dns', obj.runpath);
        end

        function readKSlices(obj, slicenums, runs)
            %READKSLICES Read in instantaneous k slices
            % Optional: numslices - only read last n slices if there are many
            
            exist("runs",'var')
            exist("numslices",'var')
            obj.kSlices = kSlice.empty;
            obj.nSlices = 0;
            if nargin < 3 || isempty(runs)
                runs = obj.run;
            end
            paths2read={};
            slicenums2read=[];
            time2read = [];
            if isempty(runs)
                ishere = true;
                switch obj.casetype
                    case 'cpu'
                        slicetime = readmatrix(fullfile(obj.runpath,'slice_time.txt'));
                        slices = dir(fullfile(obj.runpath,'kcu2_1_*'));
                    case 'gpu'
                        slicetime = readmatrix(fullfile(obj.runpath,'kslice_time.txt'));
                        slices = dir(fullfile(obj.runpath,'kcut_1_*'));
                end
                for i=1:length(slices)
                        inds(i) = str2num(slices(i).name(8:end));
                end
                [~,inds] = sort(inds);
                slices = slices(inds);
                for i=1:length(slices)
                    slicenum = str2num(slices(i).name(8:end));
                    slicenums2read(end+1) = slicenum;
                    paths2read{i} = obj.runpath;
                    for j=1:size(slicetime,1)
                        if slicetime(j,1) == slicenum
                            time2read(end+1) = slicetime(j,2);
                        end
                    end
                end

            else
                ishere=false;
                for run=runs
                    runpath = fullfile(obj.casepath,['run' num2str(run)]);
                    switch obj.casetype
                        case 'cpu'
                            slicetime = readmatrix(fullfile(runpath,'slice_time.txt'));
                            slices = dir(fullfile(runpath,'k_cuts','kcu2_1_*'));
                        case 'gpu'
                            slicetime = readmatrix(fullfile(runpath,'kslice_time.txt'));
                            slices = dir(fullfile(runpath,'k_cuts','kcut_1_*'));
                    end
                    
                    %obj.nSlices = obj.nSlices + length(slices);
                    inds = [];
                    for i=1:length(slices)
                        inds(i) = str2num(slices(i).name(8:end));
                    end
                    [~,inds] = sort(inds);
                    slices = slices(inds);
                    for i=1:length(slices)
                        paths2read{end+1} = runpath;
                        slicenum = str2num(slices(i).name(8:end));
                        slicenums2read(end+1) = slicenum;
                        for j=1:size(slicetime,1)
                            if slicetime(j,1) == slicenum
                                time2read(end+1) = slicetime(j,2);
                            end
                        end
                    end
                end
            end

            if exist("slicenums",'var')
                if length(slicenums) == 1
                    slicenums2read = slicenums2read(end-slicenums+1:end);
                    paths2read = paths2read(end-slicenums+1:end);
                    time2read = time2read(end-slicenums+1:end);
                else
                    slicenums2read = slicenums2read(slicenums);
                    paths2read = paths2read(slicenums);
                    time2read = time2read(slicenums);
                end
            end

            obj.nSlices = length(slicenums2read);

            for i=1:obj.nSlices
                fprintf('Reading slice %d/%d\n',[i obj.nSlices])
                obj.kSlices(i) = kSlice(obj.blk,obj.gas,obj.bcs,paths2read{i},slicenums2read(i),time2read(i),obj.casetype,ishere);
                %slices(i).time = slicetime(i,2);
            end
            
        end

        function readJSlices(obj, slicenums, runs)
            %READKSLICES Read in instantaneous j slices
            % Optional: numslices - only read last n slices if there are many
            
            exist("runs",'var')
            exist("numslices",'var')
            obj.jSlices = jSlice.empty;
            obj.nSlices = 0;
            if nargin < 3 || isempty(runs)
                runs = obj.run;
            end
            
            [paths2read, slicenums2read, time2read, ishere] = obj.getSlicePaths(runs);

            if exist("slicenums",'var')
                if length(slicenums) == 1
                    slicenums2read = slicenums2read(end-slicenums+1:end);
                    paths2read = paths2read(end-slicenums+1:end);
                    time2read = time2read(end-slicenums+1:end);
                else
                    slicenums2read = slicenums2read(slicenums);
                    paths2read = paths2read(slicenums);
                    time2read = time2read(slicenums);
                end
            end

            obj.nSlices = length(slicenums2read);

            for i=1:obj.nSlices
                fprintf('Reading slice %d/%d\n',[i obj.nSlices])
                obj.jSlices(i) = jSlice(obj.blk,obj.gas,paths2read{i},slicenums2read(i),time2read(i),obj.casetype,ishere);
                %slices(i).time = slicetime(i,2);
            end
            
        end

        function [paths2read, slicenums2read, time2read, ishere] = getSlicePaths(obj, runs)

            if nargin < 2 || isempty(runs)
                runs = [];
            end

            paths2read={};
            slicenums2read=[];
            time2read = [];

            if isempty(runs)
                ishere = true;
                switch obj.casetype
                    case 'cpu'
                        slicetime = readmatrix(fullfile(obj.runpath,'slice_time.txt'));
                        slices = dir(fullfile(obj.runpath, ['kcu2_1_*']));
                    case 'gpu'
                        slicetime = readmatrix(fullfile(obj.runpath,'kslice_time.txt'));
                        slices = dir(fullfile(obj.runpath,'kcut_1_*'));
                end
                for i=1:length(slices)
                        inds(i) = str2num(slices(i).name(8:end));
                end
                [~,inds] = sort(inds);
                slices = slices(inds);
                for i=1:length(slices)
                    slicenum = str2num(slices(i).name(8:end));
                    slicenums2read(end+1) = slicenum;
                    paths2read{i} = obj.runpath;
                    for j=1:size(slicetime,1)
                        if slicetime(j,1) == slicenum
                            time2read(end+1) = slicetime(j,2);
                        end
                    end
                end

            else
                ishere=false;
                for run=runs
                    runpath = fullfile(obj.casepath,['run' num2str(run)]);
                    switch obj.casetype
                        case 'cpu'
                            slicetime = readmatrix(fullfile(runpath,'slice_time.txt'));
                            slices = dir(fullfile(runpath,'k_cuts','kcu2_1_*'));
                        case 'gpu'
                            slicetime = readmatrix(fullfile(runpath,'kslice_time.txt'));
                            slices = dir(fullfile(runpath,'k_cuts','kcut_1_*'));
                    end
                    
                    %obj.nSlices = obj.nSlices + length(slices);
                    inds = [];
                    for i=1:length(slices)
                        inds(i) = str2num(slices(i).name(8:end));
                    end
                    [~,inds] = sort(inds);
                    slices = slices(inds);
                    for i=1:length(slices)
                        paths2read{end+1} = runpath;
                        slicenum = str2num(slices(i).name(8:end));
                        slicenums2read(end+1) = slicenum;
                        for j=1:size(slicetime,1)
                            if slicetime(j,1) == slicenum
                                time2read(end+1) = slicetime(j,2);
                            end
                        end
                    end
                end
            end
        end

        function readInstFlow(obj)
            %READINSTFLOW Read in instantaneous 3D flow
            obj.instFlow = volFlow(obj.blk,obj.gas,obj.bcs,obj.runpath,obj.casetype,obj.if_rans);
            
        end

        function readStartEndFlow(obj)

            fid = fopen(fullfile(obj.runpaths{1},'mean_time.txt'),'r');
            while ~feof(fid) % Use lastest mean files
                temp=fgetl(fid);
            end
%                 temp = fgetl(fid);
            fclose(fid);
            temp = str2num(temp);
            nMeanStart = temp(1)-1

            fid = fopen(fullfile(obj.runpaths{end},'mean_time.txt'),'r');
            while ~feof(fid) % Use lastest mean files
                temp=fgetl(fid);
            end
%                 temp = fgetl(fid);
            fclose(fid);
            temp = str2num(temp);
            nMeanEnd = temp(1)

            obj.startFlow = spanAveSlice(obj.runpaths{1},obj.blk,obj.gas,obj.bcs,obj.casetype,nMeanStart);
            obj.endFlow = spanAveSlice(obj.runpaths{end},obj.blk,obj.gas,obj.bcs,obj.casetype,nMeanEnd);

        end

        function writeInstFlow(obj, path)
            %WRITEINSTFLOW Write out intantaneous 3D flow to path

            if nargin < 2
                path = obj.casepath;
            end

            obj.instFlow.writeFlow(path);
        end

        function readMeanFlow(obj)
            %READMEANFLOW Read in 2D mean flow
            if strcmp(obj.runpath, obj.casepath)
                ishere = true;
            else
                ishere = false;
            end
            try
                obj.meanFlow = meanSlice(obj.blk,obj.gas,obj.bcs,obj.runpath,obj.casetype,ishere);
            catch
                try
                    data = load(fullfile(obj.runpath, 'slicesAve.mat'), 'slicesAve');
                    obj.meanFlow = data.slicesAve;
                    fprintf('Read meanSlice computed from averaging slices\n')
                catch
                    fprintf('Could not find mean flow\n')
                end
            end

        end

        function inst2mean(obj,slice)
            obj.meanFlow = meanSlice(obj.blk,obj.gas,obj.bcs,[]);
            obj.meanFlow.ro = slice.ro;
            obj.meanFlow.u = slice.u;
            obj.meanFlow.v = slice.v;
            obj.meanFlow.w = slice.w;
            obj.meanFlow.Et = slice.Et;
            obj.meanFlow.pbar = slice.p;
            obj.meanFlow.Tbar = slice.T;
            obj.meanFlow.getBCs;
        end

        function readMeanFlows(obj, calc_s_unst)
            %READMEANFLOWS Read in and average 2D mean flows for multiple
            %runs
            if nargin < 2
                calc_s_unst = false;
            end
            regions = obj.getIntRegions;
            for ir = 1:length(obj.run)
                fprintf('Reading meanSlice %d/%d (run %d)\n', [ir length(obj.run) obj.run(ir)])
                mF = meanSlice(obj.blk,obj.gas,obj.bcs,obj.runpaths{ir},obj.casetype,false, false);
                if ir == 1
                    obj.meanFlow = mF;
                    if calc_s_unst
                        [e_unst_s, ~, ~, ~, ~] = entropy_balance(mF, regions);
                    end
                else
                    obj.meanFlow.addSlice(mF)
                end
                if ir == length(obj.run)
                    if calc_s_unst
                        [e_unst_e, ~, ~, ~, ~] = entropy_balance(mF, regions);
                    end
                end
                clear mF
            end
            if calc_s_unst
                obj.e_unst = e_unst_e - e_unst_s;
            end
        end

        function readPsMeanFlows(obj, calc_s_unst)
            %READMEANFLOWS Read in and average 2D mean flows for multiple
            %runs
            if nargin < 2
                calc_s_unst = false;
            end
            regions = obj.getIntRegions;
            for ir = 1:length(obj.run)
                fprintf('Reading meanSlice %d/%d (run %d)\n', [ir length(obj.run) obj.run(ir)])
                mF = meanSlice(obj.blk,obj.gas,obj.bcs,obj.runpaths{ir},obj.casetype,false, true);
                if ir == 1
                    obj.meanFlowPS = mF;
                    if calc_s_unst
                        [e_unst_s, ~, ~, ~, ~] = entropy_balance(mF, regions);
                    end
                else
                    obj.meanFlowPS.addSlice(mF)
                end
                if ir == length(obj.run)
                    if calc_s_unst
                        [e_unst_e, ~, ~, ~, ~] = entropy_balance(mF, regions);
                    end
                end
                clear mF
            end
            if calc_s_unst
                obj.e_unst = e_unst_e - e_unst_s;
            end
        end

        function readInflowTurb(obj)
            %READINFLOWTURB Read inflow turbulence file
            obj.inflowTurb = volTurbulence(obj.casepath,obj.y_inlet,obj.bcs.ilength,obj.blk.nk,obj.bcs.lturb,obj.solver.span);
        end

        function readProbes(obj)
            
            f = fopen(fullfile(obj.runpath,'probe.txt'), 'r');
            temp = str2num(char(split(fgetl(f))));
            obj.nProbes = temp(1);
            obj.nSkip = temp(2);
            fclose(f);

%             for nProbe = 1:obj.nProbes
%                 temp = str2num(char(split(fgetl(f))));
%                 blkData.nb = temp(1);
%                 blkData.i = temp(2);
%                 blkData.j = temp(3);
%                 blkData.k = temp(4);
%                 blkData.x = obj.blk.x{blkData.nb}(blkData.i,blkData.j);
%                 blkData.y = obj.blk.y{blkData.nb}(blkData.i,blkData.j);
%                 obj.probes{nProbe} = dnsProbe(obj.runpath, nProbe, obj.nSkip, blkData, obj.gas);
%             end

            
            [obj.probes, obj.nProbes, obj.nSkip] = obj.readRunProbes(obj.runpaths{1});
            if length(obj.run) > 1
                for ir = 2:length(obj.run)

                    [probesNow, nProbesNow, nSkipNow] = obj.readRunProbes(obj.runpaths{ir});
                    if nProbesNow ~= obj.nProbes
                        fprintf('nProbes mismatch\n');
                        return
                    elseif nSkipNow ~= obj.nSkip
                        fprintf('nSkip mismatch\n');
                        return
                    else
                        for ip = 1:obj.nProbes
                            obj.probes{ip}.concatenate(probesNow{ip});
                        end
                    end
                end
            end
        end
        
        function [probes, nProbes, nSkip] = readRunProbes(obj, path)
            path
            f = fopen(fullfile(path,'probe.txt'), 'r');
            temp = str2num(char(split(fgetl(f))));
            nProbes = temp(1);
            nSkip = temp(2);
            for nProbe = 1:obj.nProbes
                temp = str2num(char(split(fgetl(f))));
                blkData.nb = temp(1);
                blkData.i = temp(2);
                blkData.j = temp(3);
                blkData.k = temp(4);
                blkData.x = obj.blk.x{blkData.nb}(blkData.i,blkData.j);
                blkData.y = obj.blk.y{blkData.nb}(blkData.i,blkData.j);
                probes{nProbe} = dnsProbe(path, nProbe, nSkip, blkData, obj.gas);
            end
            fclose(f);
        end

        function addProbe(obj, nb, i, j, k)
            blkData.nb = nb;
            blkData.i = i;
            blkData.j = j;
            if nargin < 5
                blkData.k = obj.solver.nk/2;
            else
                blkData.k = k;
            end
            blkData.x = obj.blk.x{blkData.nb}(blkData.i,blkData.j);
            blkData.y = obj.blk.y{blkData.nb}(blkData.i,blkData.j);
            obj.probes{end+1} = dnsProbe([], obj.nProbes+1, obj.nSkip, blkData, obj.gas);
            
            if isempty(obj.nProbes)
                obj.nProbes = 1;
            else
                obj.nProbes = obj.nProbes + 1;
            end
        end

        function addProbes(obj, newProbes)
            for i = 1:length(newProbes)
                blkData.nb = newProbes{i}.nb;
                blkData.i = newProbes{i}.i;
                blkData.j = newProbes{i}.j;
                blkData.k = newProbes{i}.k;
                blkData.x = obj.blk.x{blkData.nb}(blkData.i,blkData.j);
                blkData.y = obj.blk.y{blkData.nb}(blkData.i,blkData.j);
                obj.probes{end+1} = dnsProbe([], obj.nProbes+1, obj.nSkip, blkData, obj.gas);
            end
            if isempty(obj.nProbes)
                obj.nProbes = length(newProbes);
            else
                obj.nProbes = obj.nProbes + length(newProbes);
            end
        end

        function removeProbes(obj, rmProbes)
            for i = rmProbes
                obj.probes(i) = [];
                obj.nProbes = obj.nProbes-1;
            end
        end

        function plot_probes(obj, nProbe, ax)
            if nargin<3 || isempty(ax)
                ax = gca;
            end
            
            if nargin < 2 || isempty(nProbe)
                probes_to_plot = 1:obj.nProbes;
            else
                probes_to_plot = nProbe;
            end

            hold on
            C = colororder;
            sz = 25;
            for ip = probes_to_plot
                scatter(ax, obj.probes{ip}.x, obj.probes{ip}.y, sz, C(mod(ip-1,7)+1,:), 'filled')
            end
        end

        function writeProbeInput(obj, path)
            if nargin < 2
                path = obj.casepath;
            end
            f = fopen(fullfile(path, 'probe.txt'), 'w');

            fprintf(f,'%d %d\n',[obj.nProbes obj.nSkip]);
            for ip = 1:obj.nProbes
                fprintf(f,'%d %d %d %d\n', [obj.probes{ip}.nb obj.probes{ip}.i obj.probes{ip}.j obj.probes{ip}.k]);
            end
            
            fclose(f);
        end

        function s = kPlot(obj,slice,prop, varargin) % ax,lims,label,viewarea,nrepeats, rot)


            defaultAx = gca;
            defaultLims = [];
            defaultLabel = [];
            if ~isempty(obj.blk.viewarea)
                defaultViewArea = obj.blk.viewarea;
            else
                defaultViewArea = [];
            end
            if isfield(obj.blk, "n_pitchwise_repeats")
                defaultRepeats = obj.blk.n_pitchwise_repeats;
            else
                defaultRepeats = 1;
            end
            defaultRot = 0;

            p = inputParser;

            addParameter(p, 'ax', defaultAx);
            addParameter(p, 'lims', defaultLims);
            addParameter(p, 'label', defaultLabel);
            addParameter(p, 'viewarea', defaultViewArea);
            addParameter(p, 'nrepeats', defaultRepeats);
            addParameter(p, 'rot', defaultRot);
            addParameter(p, 'Interpreter', 'none');
            addParameter(p, 'ColorMap', 'parula')

            parse(p, varargin{:});

            lims = p.Results.lims;
            label = p.Results.label;
            repeats = p.Results.nrepeats;


            R = [cosd(p.Results.rot) -sind(p.Results.rot); sind(p.Results.rot) cosd(p.Results.rot)];
            ax = p.Results.ax(1);

            f = gcf;

            if ischar(prop) || isstring(prop)
                storedq = true;
                if isempty(slice)
                    q = obj.(prop);
                else
                    if strcmp(prop, "overlay")
                        q = slice.schlieren;
                        q2 = slice.vortZ;
                    else
                        q = slice.(prop);
                    end
                end
            else
                storedq = false;
                q = prop;
            end

            
            offset = 0;
            if repeats > 2
                offset = -obj.pitch;
            end
            axes(ax);
            hold on
            if size(q) == size(slice.yBL)

                xnow = slice.xO;
                ynow = slice.yO;
                ni = size(xnow,1);

                coords = [reshape(xnow, 1, []); reshape(ynow, 1, [])];
                coords = R' * coords;

                xnow = reshape(coords(1,:), ni, []);
                ynow = reshape(coords(2,:), ni, []);


                s = pcolor(ax, xnow, ynow, q);

            
            else
            for ir = 1:repeats
            for i=1:obj.NB

                xnow = obj.blk.x{i};
                ynow = obj.blk.y{i}+offset+(ir-1)*obj.pitch;
                ni = size(xnow,1);

                coords = [reshape(xnow, 1, []); reshape(ynow, 1, [])];
                coords = R' * coords;

                xnow = reshape(coords(1,:), ni, []);
                ynow = reshape(coords(2,:), ni, []);


                s = pcolor(ax, xnow, ynow, q{i});
            end
            end
            end
            shading('interp')
            axis equal
            if ~isempty(p.Results.viewarea)
                aspect = [(p.Results.viewarea(2)-p.Results.viewarea(1)) (p.Results.viewarea(4)-p.Results.viewarea(3)) 1];
                pbaspect(aspect)
                axis(p.Results.viewarea);
            end

            if storedq
                if string(prop) == "schlieren" || string(prop) == "overlay" 
                    colormap(ax,gray)
                    map = colormap(ax);
                    map = flip(map,1);
                    colormap(ax,map);
                    if isempty(label)
                        label = '$|\nabla \rho|/\rho$';
                    end
                elseif ismember(string(prop),["vortZ","v","w", "advK"])
                    val = max(abs(caxis));
                    caxis(ax,[-val val]);
                    colormap(redblue)
                elseif string(prop) == "M" && string(p.Results.ColorMap) == "redblue"
                    if isempty(lims)
                        lims = [0 2];
                    end
                    colormap(ax,redblue)
                end
            end
               
            colorbar('off');
            cb = colorbar;
            if string(prop) == "overlay"
                if size(lims, 1) == 2
                    caxis(ax,lims(2, :));
                else
                    caxis(ax,[0 100]);
                end
            elseif ~isempty(lims)
                caxis(ax,lims);
                
            end
            
            if ~isempty(label)
                cb.Label.Interpreter = 'latex';
                cb.Label.String = label;
            end
            
            if string(prop) == "overlay" 
                axes(ax);
                axis off
                area = axis(ax);
                axis on
                
                if length(p.Results.ax) > 1
                    ax2 = p.Results.ax(2);
                else
                    h = findobj('Type', 'Axes');
                    if length(h) > 1
                        ax2 = h(h~=ax);
                        cla(ax2);
                    else
                        ax2 = axes(f);
                    end
                end

                axes(ax2);
                cla(ax2);
                hold on

                a =2000;
                b = 2000;

                maxnow = 0;

                for ir = 1:repeats
                    for i=1:obj.NB
        
                        xnow = obj.blk.x{i};
                        ynow = obj.blk.y{i}+offset+(ir-1)*obj.pitch;
                        ni = size(xnow,1);
        
                        coords = [reshape(xnow, 1, []); reshape(ynow, 1, [])];
                        coords = R' * coords;
        
                        xnow = reshape(coords(1,:), ni, []);
                        ynow = reshape(coords(2,:), ni, []);
        
                        om = abs(q2{i});
                        mask = 0.5*(1+tanh((om-a)/b));
        
                        s = pcolor(ax2, xnow, ynow, q2{i});
                        alpha(s, mask);
                        maxnow = max(maxnow, max(abs(q2{i}), [], 'all'));
                    end
                end

                shading('interp')
                axis equal
                % area = axis(ax);
                axis(area);
                aspect = [(area(2)-area(1)) (area(4)-area(3)) 1];
                pbaspect(aspect);
                axis off
                colormap(ax2, redblue);
                
                if ~isempty(lims)
                    caxis(ax2, lims(1,:));
                else
                    caxis(ax2, [-maxnow maxnow]);
                end


                set(ax2, 'Position', get(ax, 'Position'));

                f.SizeChangedFcn = @copy_fig_size;

            end

            function copy_fig_size(src, ~)
                a2 = src.Children(end);
                a1 = src.Children(1);
                set(a2, 'Position', get(a1, 'Position'));
            end


        end

        function p = plotInletProf(obj,slice,prop,ax)
            if nargin < 4 || isempty(ax)
                ax = gca;
            end
            propnow = [];
            ynow = [];
            for nb = obj.blk.inlet_blocks{1}
                ynow = [ynow; obj.blk.y{nb}(1,:)];
                propnow = [propnow; slice.(prop){nb}(1,:)];
            end
            [ynow, inds] = unique(ynow);
            propnow = propnow(inds);
            p = plot(ax,propnow,ynow);

        end

        function [q,ynow] = getInletProf(obj,slice,prop)
            propnow = [];
            ynow = [];
            for nb = obj.blk.inlet_blocks{1}
                ynow = [ynow; obj.blk.y{nb}(1,:)];
                propnow = [propnow; slice.(prop){nb}(1,:)];
            end
            [ynow, inds] = unique(ynow);
            propnow = propnow(inds);
            [ynow, inds] = sort(ynow);
            propnow = propnow(inds);
            q = propnow;

        end

        function flipbook(obj,slices, prop, varargin)

            p = inputParser;

            defaultAx = gca;

            switch prop
                case 'M'
                    defaultLims = [0 1.6];
                case 'tau_w'
                    defaultLims = [0 800];
                case 'vortZ'
                    defaultLims = 1e5*[-0.7 0.7];
                case 'w'
                    defaultLims = [-3 3];
                case 'schlieren'
                    defaultLims = [0 100];
                case 'overlay'
                    defaultLims = [-0.7e5 0.7e5; 0 100];
                otherwise
                    defaultLims = [];
            end

            switch prop
                case 'M'
                    defaultLabel = '$M$';
                case 'tau_w'
                    defaultLabel = '$\tau_w$';
                case 'vortZ'
                    defaultLabel = '$\omega_z$';
                otherwise
                    defaultLabel = [];
            end

            if ~isempty(obj.blk.viewarea)
                defaultViewArea = obj.blk.viewarea;
            else
                defaultViewArea = [];
            end
            if isfield(obj.blk, "n_pitchwise_repeats")
                defaultRepeats = obj.blk.n_pitchwise_repeats;
            else
                defaultRepeats = 1;
            end
            defaultRot = 0;

            p = inputParser;

            addParameter(p, 'ax', defaultAx);
            addParameter(p, 'lims', defaultLims);
            addParameter(p, 'label', defaultLabel);
            addParameter(p, 'viewarea', defaultViewArea);
            addParameter(p, 'nrepeats', defaultRepeats);
            addParameter(p, 'rot', defaultRot);
            addParameter(p, 'Interpreter', 'none');
            addParameter(p, 'ColorMap', 'parula')

            parse(p, varargin{:});

            lims = p.Results.lims;
            label = p.Results.label;
            repeats = p.Results.nrepeats;
            ax = p.Results.ax;
            colormap = p.Results.ColorMap;

            axes(ax)
            f = gcf;
            f.UserData.slices = slices;
            f.UserData.ax = ax;
            f.UserData.homeArea = p.Results.viewarea;
            f.KeyPressFcn = @update_plot;
            f.UserData.i = 1;
            cla(ax);

            switch class(slices(1))
                case 'kSlice'
                    f.UserData.plotfcn = @(slice, area) obj.kPlot(slice,prop,'ax',ax,'lims',lims,'label',label, 'viewarea', area, 'nrepeats', repeats, 'ColorMap', colormap);
                case 'jSlice'
                    f.UserData.plotfcn = @(slice, area) obj.jPlot(slice,prop,ax,lims,label);
            end

            f.UserData.plotfcn(slices(1), p.Results.viewarea);
            title(ax, sprintf('Slice %d (%d/%d)', slices(1).nSlice, 1, length(slices)))

            function update_plot(src, event)

                area = axis(src.UserData.ax);
                switch event.Key
                    case 'leftarrow'
                        src.UserData.i = max(1, src.UserData.i - 1);
                    case 'rightarrow'
                        src.UserData.i = min(length(src.UserData.slices), src.UserData.i + 1);
                    case 'uparrow'
                        area = src.UserData.homeArea;
                end

                % cla(src.UserData.ax)
                src.UserData.plotfcn(src.UserData.slices(src.UserData.i), area)
                delete(src.UserData.ax.Children(end-(p.Results.nrepeats*obj.NB)+1:end));
                title(src.UserData.ax, sprintf('Slice %d (%d/%d)', src.UserData.slices(src.UserData.i).nSlice, src.UserData.i, length(src.UserData.slices)));

            end
        end
                

        function s = kContour(obj,slice,prop,varargin)%,ax,levels,label,fmt,linew)

            defaultAx = gca;
            defaultLevels = 1;
            defaultLabel = [];
            defaultLineW = 1;
            if ~isempty(obj.blk.viewarea)
                defaultViewArea = obj.blk.viewarea;
            else
                defaultViewArea = [];
            end
            if isfield(obj.blk, "n_pitchwise_repeats")
                defaultRepeats = obj.blk.n_pitchwise_repeats;
            else
                defaultRepeats = 1;
            end
            defaultRot = 0;

            p = inputParser;

            addParameter(p, 'ax', defaultAx);
            addParameter(p, 'levels', defaultLevels);
            addParameter(p, 'label', defaultLabel);
            addParameter(p, 'fmt', 'k');
            addParameter(p, 'LineWidth', defaultLineW);
            addParameter(p, 'FontSize', 12);
            addParameter(p, 'viewarea', defaultViewArea);
            addParameter(p, 'nrepeats', defaultRepeats);
            addParameter(p, 'rot', defaultRot);

            parse(p, varargin{:});

            levels = p.Results.levels;
            label = p.Results.label;
            repeats = p.Results.nrepeats;
            fmt = p.Results.fmt;
            linew = p.Results.LineWidth;

            R = [cosd(p.Results.rot) -sind(p.Results.rot); sind(p.Results.rot) cosd(p.Results.rot)];


            if isempty(slice)
                disp('here')
                q = obj.(prop);
            else
                q = slice.(prop);
            end

            offset = 0;
            if repeats > 2
                offset = -obj.pitch;
            end

            hold on

            if size(q) == size(slice.yBL)

                xnow = slice.xO;
                ynow = slice.yO;
                ni = size(xnow,1);

                coords = [reshape(xnow, 1, []); reshape(ynow, 1, [])];
                coords = R' * coords;

                xnow = reshape(coords(1,:), ni, []);
                ynow = reshape(coords(2,:), ni, []);


                a = smoothdata(q,1);
                a = smoothdata(a,2);
                [~,s] = contour(p.Results.ax, xnow, ynow, a, levels, fmt, 'LineWidth', linew);

            
            else
            for ir=1:repeats
            for i=1:obj.NB

                xnow = obj.blk.x{i};
                ynow = obj.blk.y{i}+offset+(ir-1)*obj.pitch;
                ni = size(xnow,1);

                coords = [reshape(xnow, 1, []); reshape(ynow, 1, [])];
                coords = R' * coords;

                xnow = reshape(coords(1,:), ni, []);
                ynow = reshape(coords(2,:), ni, []);

                a = smoothdata(q{i},1);
                a = smoothdata(a,2);
                [~,s] = contour(p.Results.ax, xnow, ynow, a, levels, fmt, 'LineWidth', linew);
            end
            end
            end
            if ~isempty(p.Results.viewarea)
                aspect = [(p.Results.viewarea(2)-p.Results.viewarea(1)) (p.Results.viewarea(4)-p.Results.viewarea(3)) 1];
                pbaspect(aspect)
                axis(p.Results.viewarea);
            end
            if string(prop) == "schlieren"
                if nargin < 6
                    %label = '$|\nabla \rho|/\rho$';
                end
            end
            if ~isempty(label)
                label;
                cb.Label.Interpreter = 'latex';
                cb.Label.String = label;
            end
            
            set(p.Results.ax, 'FontSize', p.Results.FontSize)
        end

        
        function s = BLkPlot(obj,slice,prop,ax,lims,label,viewarea,nrepeats)

            repeats = 1;
            if isfield(obj.blk, "n_pitchwise_repeats")
                repeats = obj.blk.n_pitchwise_repeats;
            end
            if nargin > 7 && ~isempty(repeats)
                repeats = nrepeats;
            end
            
            if nargin < 4 || isempty(ax)
                ax = gca;
            end

            %slice
            if isempty(slice)
                disp('here')
                q = obj.(prop);
            else
                q = slice.(prop);
            end
            hold on
            offset = 0;
            if repeats > 2
                offset = -obj.pitch;
            end
            xLE = slice.xO(1,1);
            yLE = slice.yO(1,1);
            xTE = slice.xO(end,1);
            yTE = slice.yO(end,1);

            stag = atan((yTE-yLE)/(xTE-xLE));
            R = [cos(stag) (-sin(stag)); sin(stag) cos(stag)];
            R = inv(R);
            for ir = 1:repeats
                for ib=1:slice.NB
                    xnow = [];
                    ynow = [];
                    [ni, nj] = size(obj.blk.x{ib});
                    for i=1:ni
                        for j=1:nj
                            r = R*[obj.blk.x{ib}(i,j); obj.blk.y{ib}(i,j)+(ir-1)*obj.pitch];
                            xnow(i,j) = r(1);
                            ynow(i,j) = r(2);
                        end
                    end
                    s = pcolor(ax, xnow, ynow, q{ib});
                end
            end
            shading('interp')
            axis equal
            if nargin > 6 && ~isempty(viewarea)
                aspect = [(viewarea(2)-viewarea(1)) (viewarea(4)-viewarea(3)) 1];
                pbaspect(aspect)
                axis(viewarea);
            elseif ~isempty(obj.blk.viewarea)
                pbaspect(obj.blk.aspect);
                axis(obj.blk.viewarea);
            end
            if string(prop) == "schlieren"
                colormap(gray)
                map = colormap;
                map = flip(map,1);
                colormap(map);
                if nargin < 6
                    label = '$|\nabla \rho|/\rho$';
                end
            elseif ismember(string(prop),["vortZ","v","w"])
                val = max(abs(caxis));
                caxis([-val val]);
                colormap(redblue)
            end
               
            
            if nargin > 4 && ~isempty(lims)
                caxis(lims)
            end
            if nargin > 5 || exist("label",'var')
                label;
                cb = colorbar('southoutside');
                cb.Label.Interpreter = 'latex';
                cb.Label.String = label;
            end
            
            %set(ax, 'FontSize', 12)
        end

%         function BLkPlot(obj,slice,prop,ax,lims,label)
%             if nargin < 4 || isempty(ax)
%                 ax = gca;
%             end
%             if isempty(slice)
%                 disp('here')
%                 q = obj.(prop);
%             else
%                 q = slice.(prop);
%             end
%             hold on
%             
%             xLE = slice.xO(1,1);
%             yLE = slice.yO(1,1);
%             xTE = slice.xO(end,1);
%             yTE = slice.yO(end,1);
% 
%             stag = atan((yTE-yLE)/(xTE-xLE));
%             R = [cos(stag) (-sin(stag)); sin(stag) cos(stag)];
%             R = inv(R);
%             for ib=1:slice.NB
%             xnow = [];
%             ynow = [];
%                 [ni nj] = size(obj.blk.x{ib});
%                 for i=1:ni
%                     for j=1:nj
%                         r = R*[obj.blk.x{ib}(i,j); obj.blk.y{ib}(i,j)];
%                         xnow(i,j) = r(1);
%                         ynow(i,j) = r(2);
%                     end
%                 end
%                 pcolor(ax, xnow, ynow, q{ib});
%             end
%             shading('interp')
%             
%             axis equal
%             pbaspect([5 2 1])
% 
%             %axis([0.4 0.9 0 0.2])
%             cb = colorbar('southoutside');
%             if nargin > 4 && ~isempty(lims)
%                 caxis(lims)
%             end
%             nargin
%             if nargin > 5 && ~isempty(label)
%                 label;
%                 %cb.Label.Interpreter = 'latex';
%                 cb.Label.String = label;
%             end
% 
%             set(ax, 'FontSize', 12)
%         end

        function jPlot(obj,slice,prop,ax,lims,label)
            
            if nargin < 4 || isempty(ax)
                ax = gca;
            end
            if nargin < 5 || isempty(lims)
                lims = [];
            end
            if nargin < 6 || isempty(label)
                label = [];
            end
            if string(prop) == "cf" && isempty(slice.pdyn)
                    obj.setPdyn(slice);
            end
            slice.jPlot(prop,ax,lims,label)
            
%             if isempty(slice)
%                 disp('here')
%                 q = obj.(prop);
%             else
%                 q = slice.(prop);
%             end
% 
%             pcolor(ax, slice.X, slice.Z, q);
%             
%             shading('interp')
%             %axis([-0.2 2 -0.5 0.5])
%             axis equal
%             cb = colorbar('southoutside');
%             if nargin > 4
%                 caxis(lims)
%             end
%             nargin
%             if nargin > 5
%                 label
%                 cb.Label.String = label;
%             end
        end

        function jkPlot(obj, slice, jProp, kProp, jLims, kLims)
            subplot(3,5,1:10)
            if nargin > 4
                obj.kPlot(obj.kSlices(slice),kProp,kLims)
                subplot(3,5,11:14)
                dx = abs(obj.jSlices(slice).X(end,1)-obj.jSlices(slice).X(1,1));
                dy = abs(obj.jSlices(slice).Z(1,end)-obj.jSlices(slice).Z(1,1));
                obj.jPlot(obj.jSlices(slice),jProp,jLims);
                pbaspect(gca, [dx dy dx]);
            else
                obj.kPlot(obj.kSlices(slice),kProp)
                subplot(3,5,11:14)
                dx = abs(obj.jSlices(slice).X(end,1)-obj.jSlices(slice).X(1,1));
                dy = abs(obj.jSlices(slice).Z(1,end)-obj.jSlices(slice).Z(1,1));
                obj.jPlot(obj.jSlices(slice),jProp);
                pbaspect(gca, [dx dy dx]);
            end
        end

        function p = plot_blade(obj,fmt, rot, nclock)
            % if nargin < 2 || isempty(fmt)
            %     fmt = 'k';
            % end
            if nargin < 3 || isempty(rot)
                rot = 0;
            end
            if nargin<2 || isempty(fmt)
                fmt = '';
            end

            if nargin < 4
                nclock = 0;
            end
            
            [xsurf, ysurf] = obj.get_surf;
            ysurf = ysurf + nclock*obj.blk.pitch;
            R = [cosd(rot) -sind(rot); sind(rot) cosd(rot)];
            coords = [reshape(xsurf, 1, []); reshape(ysurf, 1, [])];
            coords = R' * coords;

            if ~isempty(fmt)
            p = plot(coords(1,:),coords(2,:),fmt);
            else
            p = patch(coords(1,:), coords(2,:), [0.8 0.8 0.8]);
            p.LineWidth = 1.5;
            end
        end

        function [xsurf, ysurf] = get_surf(obj)

            xsurf = [];
            ysurf = [];
            for i=1:length(obj.blk.oblocks)
                xtemp = obj.blk.x{obj.blk.oblocks(i)}(:,end);
                ytemp = obj.blk.y{obj.blk.oblocks(i)}(:,end);
                if obj.blk.oblocks_flip(i) == 1
                    xtemp = flip(xtemp);
                    ytemp = flip(ytemp);
                end
                xsurf = [xsurf xtemp'];
                ysurf = [ysurf ytemp'];
            end

            [~, iLE] = min(xsurf);

            if ysurf(iLE) > ysurf(iLE+1)
                xsurf = [xsurf(iLE:-1:1) xsurf(end:-1:iLE+1)];
                ysurf = [ysurf(iLE:-1:1) ysurf(end:-1:iLE+1)];
            else
                xsurf = [xsurf(iLE:end) xsurf(1:iLE-1)];
                ysurf = [ysurf(iLE:end) ysurf(1:iLE-1)];
            end

        end

        function plot_surf_prop(obj,slice,prop,ax,lims)

            if nargin < 4 || isempty(ax)
                ax = gca;
            end

            q = slice.(prop);
            xsurf = [];
            for i=1:length(obj.blk.oblocks)
                temp = obj.blk.x{obj.blk.oblocks(i)}(:,end);
                qtemp = q{obj.blk.oblocks(i)}(:,end);
                if obj.blk.oblocks_flip(i) == 1
                    temp = flip(temp);
                    qtemp = flip(qtemp);
                end
                xsurf = [xsurf temp'];
                qsurf = [qsurf]
            end

            [~,iLE] = min(xsurf);
            [~,iTE] = max(xsurf);
            xsurf = xsurf(iLE:iTE);
            q = q(iLE:iTE);

            plot(ax,xsurf,q)

            if nargin == 5
                ylim(lims)
            end
        end

        function plot_BL_profile(obj,slice,x,prop,ax)
            
        end

        function setTrip(obj)
            mode = input('Trip mode:');
            obj.trip.mode = mode;
            x = input('x location:');
            switch mode
                case 1
                    obj.trip.amp = input('trip_amp:');
                    obj.trip.scale = input('trip_scale:');
                case {2, 3}
                    obj.trip.scale = input('trip_scale:');
                case 4
                    obj.trip.scale = input('ttrip_scale:');
                    obj.trip.nbumps = input('nbumps:');
            end

            points = obj.x2point(x);
            obj.trip.x = obj.blk.x{points{1}.nb}(points{1}.i,end);
            obj.trip.y = obj.blk.y{points{1}.nb}(points{1}.i,end);
            obj.trip.x2 = obj.blk.x{points{2}.nb}(points{2}.i,end);
            obj.trip.y2 = obj.blk.y{points{2}.nb}(points{2}.i,end);
            obj.iTrip = true;

        end

        function getOffsetTripCoords(obj)
            mode = input('Trip mode:');
            x = input('x location:');
            offset = input('Surface offset:');
            switch mode
                case 1
                    obj.trip.amp = input('trip_amp:');
                    obj.trip.scale = input('trip_scale:');
                case {2, 3}
                    obj.trip.scale = input('trip_scale:');
                case 4
                    obj.trip.scale = input('ttrip_scale:');
                    obj.trip.nbumps = input('nbumps:');
            end

            points = obj.x2point(x);
            x1 = obj.blk.x{points{1}.nb}(points{1}.i,end);
            y1 = obj.blk.y{points{1}.nb}(points{1}.i,end);
            x2 = obj.blk.x{points{2}.nb}(points{2}.i,end);
            y2 = obj.blk.y{points{2}.nb}(points{2}.i,end);
            
            n = obj.meanFlow.surface_normal(x1);
            x1 = x1+offset*n(1);
            y1 = y1+offset*n(2);

            disp('x1, y1:')
            fprintf('%20.16e %20.16e\n', x1, y1);
            disp('x2, y2:')
            fprintf('%20.16e %20.16e\n', x2, y2);

        end

        function point = x2point(obj, x)
            % Find closest point on top surface
            xprof = [];
            nbprof = [];
            iprof = [];
            for i=1:length(obj.blk.oblocks)
                ib = obj.blk.oblocks(i);
                if (obj.blk.next_block{ib}.jp == 0) && (obj.blk.next_patch{ib}.jp == 3)
                    xnow = obj.blk.x{ib}(:,end)';
                else
                    xnow = obj.blk.x{ib}(:,1)';
                end
                is = 1:length(xnow);
                if obj.blk.oblocks_flip(i)
                    xnow = flip(xnow);
                    is = flip(is);
                end
                xprof = [xprof xnow];
                nbprof = [nbprof ib*ones(size(xnow))];
                iprof = [iprof is];
            end
            ni = length(xprof);
            [~,iLE] = min(xprof);
            [~,iTE] = max(xprof);
            iLE = iLE+ni;
            iTE = iTE+ni;
            xprof = [xprof xprof xprof];
            iprof = [iprof iprof iprof];
            nbprof = [nbprof nbprof nbprof];
            xss = xprof(iLE:iTE);
            xps = xprof(iTE:iLE+ni);
            iss = iprof(iLE:iTE);
            ips = iprof(iTE:iLE+ni);
            nbss = nbprof(iLE:iTE);
            nbps = nbprof(iTE:iLE+ni);

            [~, i1] = min(abs(xss-x));
            [~, i2] = min(abs(xps-x));
            point{1}.i = iss(i1);
            point{1}.nb = nbss(i1);
            point{2}.i = ips(i2);
            point{2}.nb = nbps(i2);
        end

        function value = get.ftrip(obj)
            if ~obj.iTrip
                disp('Setting trip with default params')
               obj.setTrip(0.1, 1, 100)
            end
            value = cell(1,obj.NB);
            for nb = 1:obj.NB
                dist2 = ((obj.blk.x{nb} - obj.trip.x).^2 + (obj.blk.y{nb} - obj.trip.y).^2);
                if obj.trip.mode == 1
                    f = obj.trip.amp*exp(-dist2*obj.trip.scale);
                    vagenerateue{nb} = min(f,1.0);
                else
                    value{nb} = double(sqrt(dist2)<obj.trip.scale);
                end

            end
        end

        function clearTrip(obj)

            obj.trip.amp
            obj.trip.scale
            [obj.trip.nb, obj.trip.i]
            obj.trip.x
            obj.trip.y
            obj.iTrip = false;
        end

        function plotTripFactor(obj,ax)
            if nargin < 2 || isempty(ax)
                ax = gca;
            end
            if ~obj.iTrip
                disp('Must setup trip first')
            else
                ynow = obj.meanFlow.yBL(obj.meanFlow.x2ind(obj.trip.x),:);
                f = obj.ftrip{obj.trip.nb}(obj.trip.i,end:-1:1);
                U = obj.meanFlow.U(obj.meanFlow.x2ind(obj.trip.x),:);
                if nargin > 1
                    plot(ax, f, ynow)
                    hold on
                    plot(ax, U/max(U), ynow)
                else
                    plot(f, ynow)
                    hold on
                    plot(U/max(U), ynow)
                end
            end
        end

        function readInletProfile(obj)

            obj.inletProf = inletProfile(obj.blk, obj.gas, obj.bcs, obj.runpath);

        end

        function writeInletProfile(obj)

            obj.inletProf.writeProf(obj.runpath);

        end

        function readTrip(obj)
            filepath = fullfile(obj.casepath,'trip.txt');
            if ~exist(filepath)
                disp("Can't find trip.txt")
            else
                fid = fopen(filepath);
                obj.trip.mode = str2num(fgetl(fid));
                temp = fgetl(fid);
                temp = str2num(temp);
                switch obj.trip.mode
                    case 1
                        obj.trip.amp = temp(1);
                        obj.trip.scale = temp(2);
                        xtmp = temp(3);
                    case {2,3}
                        obj.trip.scale = temp(1);
                        xtmp = temp(2);
                    case {4,5}
                        obj.trip.scale = temp(1);
                        xtmp = temp(2);
                        if length(temp) < 4
                            obj.trip.nbumps = str2num(fgetl(fid));
                        else
                            obj.trip.nbumps = temp(4);
                        end
                end

                fclose(fid);
                point = obj.x2point(xtmp);
                obj.trip.nb = point{1}.nb; obj.trip.i = point{1}.i;
                obj.trip.x = obj.blk.x{obj.trip.nb}(obj.trip.i,end);
                obj.trip.y = obj.blk.y{obj.trip.nb}(obj.trip.i,end);
                obj.iTrip = true;
            end
        end

        function writeTecplotFiles(obj)
            obj.instFlow.write_tecplot_files(obj.runpath);
        end

        function writeTripInput(obj, dir)
            if ~obj.iTrip
                disp('Must setup trip first')
            else
                if nargin > 1
                    fid = fopen(fullfile(dir, 'trip.txt'), 'w');
                else
                    fid = fopen(fullfile(obj.casepath, 'trip.txt'), 'w');
                end
                switch obj.trip.mode
                    case 1
                        fprintf(fid, '%d\n', obj.trip.mode);
                        fprintf(fid,'%20.16e %20.16e %20.16e %20.16e\n', obj.trip.amp, obj.trip.scale, obj.trip.x, obj.trip.y);
                    case {2, 3}
                        fprintf(fid, '%d\n', obj.trip.mode);
                        fprintf(fid,'%20.16e %20.16e %20.16e %20.16e\n', obj.trip.scale, obj.trip.x, obj.trip.y);
                    case 4
                        fprintf(fid, '%d\n', obj.trip.mode);
                        fprintf(fid,'%20.16e %20.16e %20.16e %d\n', obj.trip.scale, obj.trip.x, obj.trip.y, obj.trip.nbumps);
                    case 5
                        fprintf(fid, '%d\n', obj.trip.mode);
                        fprintf(fid,'%20.16e %20.16e %20.16e %d\n', obj.trip.scale, obj.trip.x, obj.trip.y, obj.trip.nbumps);
                        fprintf(fid,'%20.16e %20.16e\n', obj.trip.x2, obj.trip.y2);
                end
                fclose(fid);
            end
        end

        function readRANSSlice(obj,turb,mod)
            ransdir = fullfile('RANS','cwl90');
            obj.RANSSlices{turb.mod} = RANSSlice(ransdir,data,obj.blk,obj.gas);
        end

        function [e_unst, e_conv, e_diss, e_irrev, e_rev, e_n] = entropy_budget(obj, normalise)

            if nargin < 2
                normalise = false;
            end

            regions = obj.getIntRegions;
            e_unst = (obj.area_integral(obj.endFlow.ros, regions) - ...
                obj.area_integral(obj.startFlow.ros, regions))/obj.meanFlow.meanTime;

            for ib = 1:obj.NB
                [DrousDx, ~] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.meanFlow.rous{ib});
                [~, DrovsDy] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.meanFlow.rovs{ib});
                conv_prop{ib} = DrousDx + DrovsDy;

                [DqxDx,~] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.meanFlow.rev_gen_x{ib});
                [~, DqyDy] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.meanFlow.rev_gen_y{ib});
                rev_prop{ib} = -(DqxDx + DqyDy);
            end

            e_conv = obj.area_integral(conv_prop, regions);
            e_diss = obj.area_integral(obj.meanFlow.diss_T, regions);
            e_irrev = obj.area_integral(obj.meanFlow.irrev_gen, regions);
            e_rev = obj.area_integral(rev_prop, regions);
            e_n = (e_conv + e_unst) - (e_diss + e_irrev + e_rev);

            if normalise
                for i=1:length(regions)
                    e_unst(i) = e_unst(i)/abs(e_conv(i));
                    e_conv(i) = e_conv(i)/abs(e_conv(i));
                    e_diss(i) = e_diss(i)/abs(e_conv(i));
                    e_irrev(i) = e_irrev(i)/abs(e_conv(i));
                    e_rev(i) = e_rev(i)/abs(e_conv(i));
                    e_n(i) = e_n(i)/abs(e_conv(i));
                end
            else
                factor = abs(e_diss(1));
                for i=1:length(regions)
                    e_unst(i) = e_unst(i)/factor;
                    e_conv(i) = e_conv(i)/factor;
                    e_diss(i) = e_diss(i)/factor;
                    e_irrev(i) = e_irrev(i)/factor;
                    e_rev(i) = e_rev(i)/factor;
                    e_n(i) = e_n(i)/factor;
                end
            end

            groupLabels = {};
            for ir = 1:length(regions)
                stackData(ir, 1, 1) = e_conv(ir);
                stackData(ir, 1, 2) = e_unst(ir);
                stackData(ir, 1, 3:5) = 0;
            
                stackData(ir, 2, 1:2) = 0;
                stackData(ir, 2, 3) = e_diss(ir);
                stackData(ir, 2, 4) = e_irrev(ir);
                stackData(ir, 2, 5) = e_rev(ir);
                groupLabels{ir} = regions{ir}.label;
            end
            
            h = plotBarStackGroups(stackData, groupLabels);
            c = colororder;
            c = c(1:5,:);
            c = repelem(c,size(h,1),1); 
            c = mat2cell(c,ones(size(c,1),1),3);
            set(h,{'FaceColor'},c);
            legend(h(1,:),'\epsilon_S','\epsilon_{unst}','\epsilon_\phi','\epsilon_{irrev}','\epsilon_{rev}','Location','northeast')
            set(gca,'FontSize',12)
            if normalise
                ylabel('\epsilon/|\epsilon_\Phi|')
            else
                ylabel('\epsilon/|\epsilon_{\Phi, Pre-Shock}|')
            end
            set(gca, 'XTickLabelRotation',20)
        end

%         function [e_s, e_phi, e_irrev, e_N] = entropy_balance(obj)
% 
%             e_s = 0.0;
%             for ii=1:length(obj.blk.io_surfaces.blks)
%                 
%                 prop(:,:,1) = obj.meanFlow.ro{obj.blk.io_surfaces.blks(ii)}.* ...
%                     obj.meanFlow.s{obj.blk.io_surfaces.blks(ii)}.* ...
%                     obj.meanFlow.u{obj.blk.io_surfaces.blks(ii)};
% 
%                 prop(:,:,2) = obj.meanFlow.ro{obj.blk.io_surfaces.blks(ii)}.* ...
%                     obj.meanFlow.s{obj.blk.io_surfaces.blks(ii)}.* ...
%                     obj.meanFlow.v{obj.blk.io_surfaces.blks(ii)};
%                 
%                 temp(ii) = obj.surface_integral(prop, obj.blk.io_surfaces.blks(ii), ...
%                     obj.blk.io_surfaces.types(ii));
% 
%                 e_s = e_s+temp(ii);
% 
%                 if obj.blk.io_surfaces.blks(ii) == 3
%                     temp(ii)
%                 end
% 
%                 clear prop
% 
%             end
% 
%             e_phi = 0.0;
%             e_irrev = 0.0;
%             
%             parfor ib = 1:obj.NB
%                 ib
%                 prop_phi = obj.meanFlow.diss_T{ib};
%                 prop_phi = obj.meanFlow.diss{ib}./obj.meanFlow.T{ib};
%                 prop_irrev = obj.meanFlow.irrev_gen{ib};
%                 e_phi = e_phi + obj.area_integral(prop_phi, ib);
%                 e_irrev = e_irrev + obj.area_integral(prop_irrev, ib);
%             end
% 
%             e_N = e_s - e_phi - e_irrev;
% 
%             
%         end

        function value = surface_integral(obj, prop, ib, type)
            if ib == 3
                ib
                type
            end
            switch type
                case 1
                    xnow = obj.blk.x{ib}(1,end:-1:1);
                    ynow = obj.blk.y{ib}(1,end:-1:1);
                    qx = prop(1,end:-1:1,1);
                    qy = prop(1,end:-1:1,2);
                case 2 
                    xnow = obj.blk.x{ib}(end,:);
                    ynow = obj.blk.y{ib}(end,:);
                    qx = prop(end,:,1);
                    qy = prop(end,:,2);
                case 3
                    xnow = obj.blk.x{ib}(:,1);
                    ynow = obj.blk.y{ib}(:,1);
                    qx = prop(:,1,1);
                    qy = prop(:,1,2);
                case 4
                    xnow = obj.blk.x{ib}(end:-1:1,end);
                    ynow = obj.blk.y{ib}(end:-1:1,end);
                    qx = prop(end:-1:1,end,1);
                    qy = prop(end:-1:1,end,2);
            end

            xnow = reshape(xnow,1,[]);
            ynow = reshape(ynow,1,[]);
            qx = reshape(qx,1,[]);
            qy = reshape(qy,1,[]);

            dx = xnow(2:end)-xnow(1:end-1);
            dy = ynow(2:end)-ynow(1:end-1);
            ds = sqrt(dx.^2+dy.^2);
            s = [0 cumsum(ds)];
            
            n(1,:) = [dy(1)/ds(1) -dx(1)/ds(1)];
            n(length(xnow),:) = [dy(end)/ds(end) -dx(end)/ds(end)];

            for i=2:length(xnow)-1
                n1 = [dy(i-1)/ds(i-1) -dx(i-1)/ds(i-1)];
                n2 = [dy(i)/ds(i) -dx(i)/ds(i)];
                n(i,:) = (n1+n2)/norm(n1+n2);
            end
            

            integrand = n(:,1).*qx' + n(:,2).*qy';
            value = trapz(s',integrand);


        end

        function value = area_integral(obj, prop, regions)
            value = [];
            for r = [regions{:}]

                im = r.is:r.ie-1;
                ip = r.is+1:r.ie;
                jm = r.js:r.je-1;
                jp = r.js+1:r.je;
                
                x1 = obj.blk.x{r.nb}(im,jm);
                x2 = obj.blk.x{r.nb}(ip,jm);
                x3 = obj.blk.x{r.nb}(ip,jp);
                x4 = obj.blk.x{r.nb}(im,jp);

                y1 = obj.blk.y{r.nb}(im,jm);
                y2 = obj.blk.y{r.nb}(ip,jm);
                y3 = obj.blk.y{r.nb}(ip,jp);
                y4 = obj.blk.y{r.nb}(im,jp);

                area = 0.5*( (x1.*y2 + x2.*y3 + x3.*y4 + x4.*y1) ...
                    - (x2.*y1 + x3.*y2 + x4.*y3 + x1.*y4) );

                q = 0.25*(prop{r.nb}(im,jm) + ...
                          prop{r.nb}(ip,jm) + ...
                          prop{r.nb}(ip,jp) + ...
                          prop{r.nb}(im,jp));

                value(end+1) = sum(q.*area,"all");

            end
%             for i=1:obj.blk.blockdims(nb,1)-1
%                 for j=1:obj.blk.blockdims(nb,2)-1
%                     xnow = [obj.blk.x{nb}(i,j) obj.blk.x{nb}(i+1,j) ...
%                         obj.blk.x{nb}(i+1,j+1) obj.blk.x{nb}(i,j+1)];
%                     ynow = [obj.blk.y{nb}(i,j) obj.blk.y{nb}(i+1,j) ...
%                         obj.blk.y{nb}(i+1,j+1) obj.blk.y{nb}(i,j+1)];
% 
%                     q = 0.25*(prop(i,j)+prop(i+1,j)+prop(i+1,j+1)+prop(i,j+1));
% 
%                     area = polyarea(xnow,ynow);
%                     value = value+area*q;
%                 end
%             end
        end

%         function [e_s, e_phi, e_irrev, e_rev] = entropy_balance(obj)
%             
%             pool = gcp('nocreate');
%             if isempty(pool)
%                 parpool;
%             end
%             
%             %
%             
%             celldims = obj.blk.blockdims(:,1:2) - 1;
%             cells.x = zeros([sum(prod(celldims,2)) 1]);
%             cells.y = zeros([sum(prod(celldims,2)) 1]);
%             cells.area = zeros([sum(prod(celldims,2)) 1]);
%             cells.unst = zeros([sum(prod(celldims,2)) 1]);
%             cells.conv = zeros([sum(prod(celldims,2)) 1]);
%             cells.diss = zeros([sum(prod(celldims,2)) 1]);
%             cells.irrev = zeros([sum(prod(celldims,2)) 1]);
%             cells.rev = zeros([sum(prod(celldims,2)) 1]);
%             
%             
%             %
%             
%             parfor ib=1:obj.NB
%                 nentries = celldims(ib,1)*celldims(ib,2);
%             
%             
%                 [drousdx,~] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.meanFlow.rous{ib});
%                 [~,drovsdy] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.meanFlow.rovs{ib});
%                 [dqx_Tdx,~] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.meanFlow.rev_gen_x{ib});
%                 [~,dqy_Tdy] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.meanFlow.rev_gen_y{ib});
%             
%                 unst_prop = obj.meanFlow.ro{ib}.*dsdt;
%             
%                 conv_prop = obj.meanFlow.ro{ib}.*(obj.meanFlow.u{ib}.*dsdx + ...
%                    obj.meanFlow.v{ib}.*dsdy);
% 
%                 conv_prop = drousdx + drovsdy;
%             
%                 diss_prop = obj.meanFlow.diss_T{ib};
%                 irrev_prop = obj.meanFlow.irrev_gen{ib};
%                 rev_prop = dqx_Tdx + dqy_Tdy;
%             
%                 xtmp{ib} = zeros(nentries, 1);
%                 ytmp{ib} = zeros(nentries, 1);
%                 areatmp{ib} = zeros(nentries, 1);
%                 unsttmp{ib} = zeros(nentries,1);
%                 convtmp{ib} = zeros(nentries, 1);
%                 disstmp{ib} = zeros(nentries, 1);
%                 irrevtmp{ib} = zeros(nentries, 1);
%                 revtmp{ib} = zeros(nentries, 1);
%                 
%                 for i=1:celldims(ib,1)
%                     if mod(i, 50) == 0
%                          sprintf('Block %d, i=%d/%d', ib, i, celldims(ib,1))
%                     end
%                     for j=1:celldims(ib,2)
%                         pos = celldims(ib,2)*(i-1)+j;
%                         xnow = [obj.blk.x{ib}(i,j) obj.blk.x{ib}(i+1,j) ...
%                             obj.blk.x{ib}(i+1,j+1) obj.blk.x{ib}(i,j+1)];
%                         ynow = [obj.blk.y{ib}(i,j) obj.blk.y{ib}(i+1,j) ...
%                             obj.blk.y{ib}(i+1,j+1) obj.blk.y{ib}(i,j+1)];
%             
%                         xtmp{ib}(pos) = mean(xnow);
%                         ytmp{ib}(pos) = mean(ynow);
%             
%                         area = polyarea(xnow,ynow);
%                         areatmp{ib}(pos) = area;
%             
%                         unsttmp{ib}(pos) = 0.25*(unst_prop(i,j)+unst_prop(i+1,j)+unst_prop(i+1,j+1)+unst_prop(i,j+1))*area;
%                         convtmp{ib}(pos) = 0.25*(conv_prop(i,j)+conv_prop(i+1,j)+conv_prop(i+1,j+1)+conv_prop(i,j+1))*area;
%                         disstmp{ib}(pos) = 0.25*(diss_prop(i,j)+diss_prop(i+1,j)+diss_prop(i+1,j+1)+diss_prop(i,j+1))*area;
%                         irrevtmp{ib}(pos) = 0.25*(irrev_prop(i,j)+irrev_prop(i+1,j)+irrev_prop(i+1,j+1)+irrev_prop(i,j+1))*area;
%                         revtmp{ib}(pos) = 0.25*(rev_prop(i,j)+rev_prop(i+1,j)+rev_prop(i+1,j+1)+rev_prop(i,j+1))*area;
%             
%                     end
%                 end
%             end
%             
%             for ib = 1:obj.NB
%             
%                 offset = sum(prod(celldims(1:ib,:),2));
%                 nentries = celldims(ib,1)*celldims(ib,2);
%             
%                 cells.x(offset-nentries+1:offset) = xtmp{ib};
%                 cells.y(offset-nentries+1:offset) = ytmp{ib};
%                 cells.area(offset-nentries+1:offset) = areatmp{ib};
%                 %cells.unst(offset-nentries+1:offset) = unsttmp{ib};
%                 cells.conv(offset-nentries+1:offset) = convtmp{ib};
%                 cells.diss(offset-nentries+1:offset) = disstmp{ib};
%                 cells.irrev(offset-nentries+1:offset) = irrevtmp{ib};
%                 cells.rev(offset-nentries+1:offset) = revtmp{ib};
%             end
%             
%             %%
%             
%             indices = cells.x > 0.75 & cells.x < 1.2 & cells.y > -0.1 & cells.y < 0.1;
%             
%             %e_unst = sum(cells.unst(indices));
%             e_s = sum(cells.conv(indices));
%             e_phi = sum(cells.diss(indices));
%             e_irrev = sum(cells.irrev(indices));
%             e_rev = sum(cells.rev(indices));
% 
%         end

        function plot_surf(obj)

            for ib = obj.blk.oblocks
            
                plot3(obj.blk.x{ib}(:,end), obj.blk.y{ib}(:,end), zeros(size(obj.blk.x{ib}(:,end))), 'k');
                hold on
                plot3(obj.blk.x{ib}(:,end), obj.blk.y{ib}(:,end), obj.blk.span*ones(size(obj.blk.x{ib}(:,end))), 'k');
                axis equal
            end


        end

        function plot_surf_2d(obj)

            for ib = obj.blk.oblocks
                hold on
                plot(obj.blk.x{ib}(:,end), obj.blk.y{ib}(:,end), 'k');
                axis equal
            end


        end

        function plot_mesh(obj, skip, ax, label_blocks)
            if nargin<2 || isempty(skip)
                skip=8;
            end

            if nargin<3 || isempty(ax)
                ax = gca;
            end

            if nargin < 4 || isempty(label_blocks)
                label_blocks = false;
            end
            
            hold on
            for ib=1:obj.NB
                [ni, nj] = size(obj.blk.x{ib});
                for j=[1:skip:ni ni]
                    if (j==1) || (j==ni)
                        plot(ax, obj.blk.x{ib}(j,:),obj.blk.y{ib}(j,:),'r','LineWidth',1)
                    else
                        plot(ax, obj.blk.x{ib}(j,:),obj.blk.y{ib}(j,:),'k')
                    end
                end
                for j=[1:skip:nj nj]
                    if (j==1) || (j==nj)
                        plot(ax, obj.blk.x{ib}(:,j),obj.blk.y{ib}(:,j),'r','LineWidth',1)
                    else
                        plot(ax, obj.blk.x{ib}(:,j),obj.blk.y{ib}(:,j),'k')
                    end
                end
            end

            if label_blocks
                for ib = 1:obj.NB
                    i = floor(obj.blk.blockdims(ib,1)/2);
                    j = floor(obj.blk.blockdims(ib,2)/2);
                    text(obj.blk.x{ib}(i,j), obj.blk.y{ib}(i,j), num2str(ib),"Color",'r');
                end
            end
            
            axis equal
        end

        function setPdyn(obj,slice)
            if isempty(obj.meanFlow)
                obj.readMeanFlow;
            end
            slice.pdyn = obj.meanFlow.pdyn;
        end

        function mesh_analysis(obj, skip, ax)

            if nargin<2 || isempty(skip), skip=8; end

            mesh_analysis(obj.blk, skip)

%             if nargin < 3 || isempty(ax)
%                 ax = gca;
%             end
%             C = colororder;
%             
%             npp = obj.blk.npp;
%             
%             nij_pts = 0;
%             hold on
%             for i=1:obj.NB
%                 x = obj.blk.x{i};
%                 y = obj.blk.y{i};
%                 x = x(:,:,1);
%                 y = y(:,:,1);
%                 ni = size(x,1);
%                 nj = size(x,2);
%                 nij_pts = nij_pts + ni*nj;
%                 %nk = blk{i}.nk;
%                 for j=[1:skip:ni ni]
%                     if (j==1) || (j==ni)
%                         plot(x(j,:),y(j,:),'r','LineWidth',1)
%                     else
%                         plot(x(j,:),y(j,:),'k')
%                     end
%                 end
%                 for j=[1:skip:nj nj]
%                     if (j==1) || (j==nj)
%                         plot(x(:,j),y(:,j),'r','LineWidth',1)
%                     else
%                         plot(x(:,j),y(:,j),'k')
%                     end
%                 end
%             end

        end

        function value = get.aspect_ratio(obj)

            for i=1:obj.NB
                x = obj.blk.x{i};
                y = obj.blk.y{i};
                x = x(:,:,1);
                y = y(:,:,1);
            
                dxi = x(2:end,:)-x(1:end-1,:);
                dxi(end+1,:) = dxi(end,:);
                dxj = x(:,2:end)-x(:,1:end-1);
                dxj(:,end+1) = dxj(:,end);
                
                dyi = y(2:end,:)-y(1:end-1,:);
                dyi(end+1,:) = dyi(end,:);
                dyj = y(:,2:end)-y(:,1:end-1);
                dyj(:,end+1) = dyj(:,end);
                
                dsi = sqrt(dxi.^2+dyi.^2);
                dsj = sqrt(dxj.^2+dyj.^2);
            
                % min_spacing(i) = min(min(dsi, dsj),[],'all');
                
                value{i} = max(dsi./dsj, dsj./dsi);
                % value{i}(end+1,:) = value{i}(end,:);
                % value{i}(:,end+1) = value{i}(:,end);
            end
        end

        function writeInputFiles(obj)
            write_input_files(obj.casepath,obj.blk, ...
                obj.bcs,obj.gas,obj.solver,'topology',obj.topology,'casetype',obj.casetype);
        end

        function writeGridFiles(obj)
            write_grid(obj.casepath, obj.blk)
        end

        function get_cell_areas(obj)
            obj.cell_area = {};
            fprintf('Calculating cell areas\n');
            for ib = 1:obj.NB
                fprintf('Block %d\n',ib);
                for i=1:obj.blk.blockdims(ib,1)-1
                    if mod(i,20) == 0
                        fprintf('i = %d/%d\n',[i obj.blk.blockdims(ib,1)-1])
                    end
                    for j=1:obj.blk.blockdims(ib,2)-1
                        xnow = [obj.blk.x{ib}(i,j) obj.blk.x{ib}(i+1,j) ...
                            obj.blk.x{ib}(i+1,j+1) obj.blk.x{ib}(i,j+1)];
                        ynow = [obj.blk.y{ib}(i,j) obj.blk.y{ib}(i+1,j) ...
                            obj.blk.y{ib}(i+1,j+1) obj.blk.y{ib}(i,j+1)];

                        obj.cell_area{ib}(i,j) = abs(polyarea(xnow,ynow));
                    end
                end
            end
        end

        function [e_unst, e_conv] = mass_balance(obj)
            if isempty(obj.cell_area)
                obj.get_cell_areas;
            end
            
            prop = {};
            for ib=1:obj.NB
                [drou_dx, ~] = gradHO(obj.blk.x{ib}, obj.blk.y{ib}, obj.meanFlow.ro{ib}.*obj.meanFlow.u{ib});
                [~, drov_dy] = gradHO(obj.blk.x{ib}, obj.blk.y{ib}, obj.meanFlow.ro{ib}.*obj.meanFlow.v{ib});

                conv_prop{ib} = drou_dx + drov_dy;
                
                prop{ib}(:,:,1) = obj.meanFlow.ro{ib}.*obj.meanFlow.u{ib};
                prop{ib}(:,:,2) = obj.meanFlow.ro{ib}.*obj.meanFlow.v{ib};
                
            end

            regions = obj.getIntRegions;
            [conv_vol, conv_surf] = prop_balance(prop, obj.blk, obj.cell_area, regions);
            [conv_vol2, conv_surf2] = prop_balance(conv_prop, obj.blk, obj.cell_area, regions);

        end

        function [e_unst, e_s, e_phi, e_irrev, e_rev] = entropy_balance(obj, surface_int)
            if nargin < 2
                surface_int = true;
            end

            if isempty(obj.cell_area)
                obj.get_cell_areas;
            end

            for ib = 1:obj.NB
                s_flux{ib}(:,:,1) = obj.meanFlow.rous{ib};
                s_flux{ib}(:,:,2) = obj.meanFlow.rovs{ib};
                q_T_flux{ib}(:,:,1) = obj.meanFlow.rev_gen_x{ib};
                q_T_flux{ib}(:,:,2) = obj.meanFlow.rev_gen_y{ib};
            end

            regions = obj.getIntRegions;
            [e_s_vol, e_s_surf] = prop_balance(s_flux, obj.blk, obj.cell_area, regions);
            [e_phi, ~] = prop_balance(obj.meanFlow.diss_T, obj.blk, obj.cell_area, regions);
            [e_irrev, ~] = prop_balance(obj.meanFlow.irrev_gen, obj.blk, obj.cell_area, regions);
            [e_rev_vol, e_rev_surf] = prop_balance(q_T_flux, obj.blk, obj.cell_area, regions);
            e_unst = {};
            
            if surface_int
                e_s = e_s_surf;
                e_rev = e_rev_surf;
            else
                e_s = e_s_vol;
                e_rev = e_rev_vol;
            end 
        end

        function [e_unst, e_s, e_phi, e_irrev, e_rev] = getEntropyTerms(obj, mF)

            if nargin < 2
                mF = obj.meanFlow;
            end

            regions = obj.getIntRegions;

%             if length(obj.run) == 1
            if isempty(mF.e_s)
                [~, e_s, e_phi, e_irrev, e_rev] = entropy_balance(mF, regions);
                mF.e_s = e_s;
                mF.e_phi = e_phi;
                mF.e_irrev = e_irrev;
                mF.e_rev = e_rev;
            else
                e_s = mF.e_s;
                e_phi = mF.e_phi;
                e_irrev = mF.e_irrev;
                e_rev = mF.e_rev;
            end
            if isempty(obj.e_unst)
                e_unst = zeros(size(e_s));
            else
                e_unst = obj.e_unst;
            end
%                e_unst = zeros(size(e_unst));
%             else
%                 es_unst = {};
%                 es_s = {};
%                 es_phi = {};
%                 es_irrev = {};
%                 es_rev = {};
%                 for ir = 1:length(obj.run)
%                     [es_unst{ir}, es_s{ir}, es_phi{ir}, es_irrev{ir}, es_rev{ir}] = entropy_balance(mF, regions);
%                 end
%                 for i = 1:length(regions)
%                     e_unst(i) = es_unst{end}(i)-es_unst{1}(i);
%                     e_s(i) = 0; e_phi(i) = 0; e_irrev(i) = 0; e_rev(i) = 0;
%                     for j = 1:length(obj.run)
%                         e_s(i) = e_s(i) + es_s{j}(i);
%                         e_phi(i) = e_phi(i) + es_phi{j}(i);
%                         e_irrev(i) = e_irrev(i) + es_irrev{j}(i);
%                         e_rev(i) = e_rev(i) + es_rev{j}(i);
%                     end
%                     e_s(i) = e_s(i)/length(obj.run);
%                     e_phi(i) = e_phi(i)/length(obj.run);
%                     e_irrev(i) = e_irrev(i)/length(obj.run);
%                     e_rev(i) = e_rev(i)/length(obj.run);
%                 end
%             end
        end

        function [e_Pr, e_diss] = getTkeTerms(obj, mF)

            if nargin < 2
                mF = obj.meanFlow;
            end

            regions = obj.getIntRegions;

%             if length(obj.run) == 1
            if isempty(mF.e_Pr)
                [e_Pr, e_diss] = tke_balance(mF, regions);
                mF.e_Pr = e_Pr;
                mF.e_diss = e_diss;
            else
                e_Pr = mF.e_Pr;
                e_diss = mF.e_diss;
            end

        end

        function plot_diss_regions(obj)
%             [e_unst, e_s, e_phi, e_irrev, e_rev] = obj.getEntropyTerms;
            [~, e_s, e_phi, e_irrev, e_rev] = obj.entropy_balance;
            e_unst = obj.e_unst;
            enow = e_s/e_s(1);
            X = categorical({'Pre-shock', 'Post-shock', 'Trailing edge', 'Wake'});
            X = reordercats(X,{'Pre-shock', 'Post-shock', 'Trailing edge', 'Wake'});
            bar(X,enow)
            ylabel('\epsilon_S/\epsilon_{S, Pre-shock}')
            set(gca, 'FontSize', 12)
            set(gca, 'XTickLabelRotation', 20)
        end

        function plot_entropy_stackup(obj, normalise)
            if nargin < 2
                normalise = false;
            end

%             [e_unst, e_s, e_phi, e_irrev, e_rev] = obj.getEntropyTerms;
            [~, e_s, e_phi, e_irrev, e_rev] = obj.entropy_balance;
            if ~isempty(obj.e_unst)
                e_unst = obj.e_unst
            else
                e_unst = zeros(size(e_s));
            end
            normalise

            if normalise
                for i=1:length(e_s)
                    e_unst(i) = e_unst(i)/abs(e_s(i));
                    e_phi(i) = e_phi(i)/abs(e_s(i));
                    e_irrev(i) = e_irrev(i)/abs(e_s(i));
                    e_rev(i) = e_rev(i)/abs(e_s(i));
                    e_s(i) = e_s(i)/abs(e_s(i));
                end
            else
                factor = abs(e_phi(1));
                for i=1:length(e_s)
                    e_unst(i) = e_unst(i)/factor;
                    e_phi(i) = e_phi(i)/factor;
                    e_irrev(i) = e_irrev(i)/factor;
                    e_rev(i) = e_rev(i)/factor;
                    e_s(i) = e_s(i)/factor;
                end
            end

            for ng = 1:length(e_s)
                stackData(ng, 1, 1) = e_s(ng);
                stackData(ng, 1, 2) = e_unst(ng);
                stackData(ng, 1, 3:5) = 0;
            
                stackData(ng, 2, 1:2) = 0;
                stackData(ng, 2, 3) = e_phi(ng);
                stackData(ng, 2, 4) = e_irrev(ng);
                stackData(ng, 2, 5) = e_rev(ng);
            end
            
            groupLabels = {"Pre-shock", "Post-shock", "Trailing edge", "Wake"};
            h = plotBarStackGroups(stackData, groupLabels);
            c = colororder;
            c = c(1:5,:);
            c = repelem(c,size(h,1),1); 
            c = mat2cell(c,ones(size(c,1),1),3);
            set(h,{'FaceColor'},c);
            legend(h(1,:),'\epsilon_S','\epsilon_{unst}','\epsilon_\phi','\epsilon_{irrev}','\epsilon_{rev}','Location','northeast')
            set(gca,'FontSize',12)
            if normalise
                ylabel('\epsilon/|\epsilon_\Phi|')
            else
                ylabel('\epsilon/|\epsilon_{\Phi, Pre-Shock}|')
            end
            set(gca, 'XTickLabelRotation',20)
        end

        function value = get.trip_Re_x(obj)
            if isempty(obj.trip)
                value = [];
            else
            end
        end
        
        function value = get.Re_k(obj)
            if isempty(obj.trip)
                value = [];
            else
                tripInd = obj.meanFlow.x2ind(obj.trip.x);
                inds = obj.meanFlow.BLedgeInd;
                nunow = obj.meanFlow.nu_e;
                Unow = obj.meanFlow.U;
                Ue = Unow(tripInd,inds(tripInd));
                nu_k = nunow(tripInd);
                value = Ue*obj.trip.scale/nu_k;
            end
        end

        function regions = getIntRegions(obj, iplot)
            if nargin < 2
                iplot = false;
            end
            switch obj.topology
                case 1
                    xrangePreShock = [0.1 0.45];
                    xrangePostShock = [0.55 0.9];
        
                    regions = {};
                    regions{1}.nb = 4;
                    [~, regions{1}.is] = min(abs(obj.blk.x{4}(:,1) - xrangePreShock(1)));
                    [~, regions{1}.ie] = min(abs(obj.blk.x{4}(:,1) - xrangePreShock(2)));
                    regions{1}.js = 1;
                    regions{1}.je = 192;
                    regions{1}.label = "Pre-shock";
                    
                    regions{2}.nb = 4;
                    [~, regions{2}.is] = min(abs(obj.blk.x{4}(:,1) - xrangePostShock(1)));
                    [~, regions{2}.ie] = min(abs(obj.blk.x{4}(:,1) - xrangePostShock(2)));
                    regions{2}.js = 1;
                    regions{2}.je =  192;
                    regions{2}.label = "Post-shock";
                case 2
                    xrangePreShock = [0.16 0.59];
                    xrangePostShock = [0.70 0.98];
                    xMaxWake = 1.4;
        
                    regions = {};
                    regions{1}.nb = 6;
                    [~, regions{1}.is] = min(abs(obj.blk.x{6}(:,end) - xrangePreShock(1)));
                    [~, regions{1}.ie] = min(abs(obj.blk.x{6}(:,end) - xrangePreShock(2)));
                    regions{1}.js = 1;
                    regions{1}.je = obj.blk.blockdims(6,2);
                    
                    regions{2}.nb = 6;
                    [~, regions{2}.is] = min(abs(obj.blk.x{6}(:,end) - xrangePostShock(1)));
                    [~, regions{2}.ie] = min(abs(obj.blk.x{6}(:,end) - xrangePostShock(2)));
                    regions{2}.js = 1;
                    regions{2}.je =  obj.blk.blockdims(6,2);
                    
                    regions{3}.nb = 9;
                    regions{3}.is = 1; regions{3}.ie = size(obj.blk.x{9},1);
                    regions{3}.js = 1; regions{3}.je = size(obj.blk.x{9},2);
        
                    regions{4}.nb = 10;
                    regions{4}.is = 1;
                    [~, regions{4}.ie] = min(abs(obj.blk.x{10}(:,obj.blk.blockdims(10,2)/2) - xMaxWake));
                    regions{4}.js = 1; regions{4}.je = obj.blk.blockdims(10,2);
                case 3
                    xrangePreShock = [0.1 0.45];
                    xrangePostShock = [0.55 0.9];
        
                    regions = {};
                    regions{1}.nb = 1;
                    [~, regions{1}.is] = min(abs(obj.blk.x{1}(:,1) - xrangePreShock(1)));
                    [~, regions{1}.ie] = min(abs(obj.blk.x{1}(:,1) - xrangePreShock(2)));
                    regions{1}.js = 1;
                    regions{1}.je = 192;
                    regions{1}.label = "Pre-shock";
                    
                    regions{2}.nb = 1;
                    [~, regions{2}.is] = min(abs(obj.blk.x{1}(:,1) - xrangePostShock(1)));
                    [~, regions{2}.ie] = min(abs(obj.blk.x{1}(:,1) - xrangePostShock(2)));
                    regions{2}.js = 1;
                    regions{2}.je =  192;
                    regions{2}.label = "Post-shock";
            end
            if iplot
                ax = gca;
                obj.kPlot(obj.meanFlow, 'diss','ax',ax,'label','$\phi$');
                hold on
                p = [];
                for ir = 1:length(regions)
                    nb = regions{ir}.nb;
                    irange = regions{ir}.is:regions{ir}.ie;
                    jrange = regions{ir}.js:regions{ir}.je;
                    xline = [obj.blk.x{nb}(irange,jrange(1))' obj.blk.x{nb}(irange(end),jrange(2:end)) ...
                        obj.blk.x{nb}(irange(end-1:-1:1),jrange(end))' obj.blk.x{nb}(irange(1), jrange(end-1:-1:1))];
                    yline = [obj.blk.y{nb}(irange,jrange(1))' obj.blk.y{nb}(irange(end),jrange(2:end)) ...
                        obj.blk.y{nb}(irange(end-1:-1:1),jrange(end))' obj.blk.y{nb}(irange(1), jrange(end-1:-1:1))];
                    p(ir) = plot(ax, xline, yline, 'LineWidth', 1.5);
                    labels{ir} = regions{ir}.label;
                end
                legend(p,[labels{:}]);
            end
        end

        function interpInstFlow(obj, newcase)
            for ib = 1:obj.NB
                oldFlow = volFlowBlock(obj.runpath, ib, obj.blk, obj.gas, obj.casetype);
                newFlow = oldFlow.interpOntoNewGrid(newcase, ib);
                newFlow.writeFlow(newcase.casepath)
                clear oldFlow newFlow
            end
        end

        function interpSlice(obj, slice, newcase)
            for ib = 1:obj.NB
                oldFlow = slice.slice2flowBlock(ib);
                newFlow = oldFlow.interpOntoNewGrid(newcase, ib);
                newFlow.writeFlow(newcase.casepath,newcase.casetype);
                clear oldFlow newFlow
            end
        end

        function updateBlk(obj, newblk)
            obj.blk.x = newblk.x;
            obj.blk.y = newblk.y;
            obj.blk.nk = newblk.nk;
            obj.blk.npp = newblk.npp;
            obj.blk.z = linspace(0, obj.solver.span, obj.blk.nk);
            obj.blk.blockdims = newblk.blockdims;

            y_inlet = [];
            for i=1:length(obj.blk.inlet_blocks{1})
                y_inlet = [y_inlet obj.blk.y{obj.blk.inlet_blocks{1}(i)}(1,:,1)];
                ymidnow(i) = obj.blk.y{obj.blk.inlet_blocks{1}(i)}(1,ceil(end/2));
            end
            obj.y_inlet = sort(unique(y_inlet));
            obj.nj_inlet = length(obj.y_inlet);
            obj.inlet_width = max(obj.y_inlet) - min(obj.y_inlet);
        
        end

        function nb = find_block(obj, x, y)
            nb = [];
            for ib=1:obj.NB

                [xb, yb] = block_boundary(obj.blk, ib);
                if inpolygon(x,y,xb,yb)
                    nb = ib;
                    break
                end
            end

        end

        function [i,j,nb] = find_ij(obj, x, y)
            nb = obj.find_block(x,y);
            xnow = obj.blk.x{nb};
            ynow = obj.blk.y{nb};
            ds = sqrt((xnow-x).^2 + (ynow-y).^2);
            [~,ind] = min(ds,[],'all');
            [i,j] = ind2sub(size(xnow),ind);
        end

        function newProbes = interp_probes(obj, oldCase)
            oldProbes = oldCase.probes;
            newProbes = {};
            for ip = 1:oldCase.nProbes
                [i,j,nb] = obj.find_ij(oldProbes{ip}.x,oldProbes{ip}.y);
                blkData.i = i;
                blkData.j = j;
                blkData.nb = nb;
                blkData.k = ceil(obj.solver.nk/2);
                blkData.x = obj.blk.x{oldProbes{ip}.nb}(i,j);
                blkData.y = obj.blk.y{oldProbes{ip}.nb}(i,j);
                newProbes{ip} = dnsProbe([],ip,obj.nSkip,blkData,obj.gas);
            end
        end

        function newCase = instantiate(obj)
            newCase = DNS_case;
        end

        function newCase = setup_new_case(obj, name, newBlk, iWrite)
            if nargin < 4
                iWrite = false;
            end
            newCase = obj.instantiate;
            newCase.casename=name;
            newCase.casetype = obj.casetype;
            newCase.topology = obj.topology;
            newCase.casepath = fullfile(fileparts(obj.casepath), name);
            if ~exist(newCase.casepath,'dir')
                mkdir(newCase.casepath);
            end
            newCase.blk = newBlk;
            fields = {'next_block', 'next_patch', 'corner'};
            for i = 1:length(fields)
                if ~isfield(newBlk, fields{i})
                    newCase.blk.(fields{i}) = obj.blk.(fields{i});
                end
            end
            newCase.solver = obj.solver;
            newCase.solver.nk = newBlk.nk;
            newCase.bcs = obj.bcs;
            newCase.gas = obj.gas;
            newCase.NB = length(newBlk.x);
            newCase.trip = obj.trip;
            newCase.iTrip = obj.iTrip;
            newCase.construct_o_blocks;
            newCase.compute_blk_metadata;
            newCase.casetype = obj.casetype;

            if iWrite
                newCase.writeInputFiles;
                newCase.writeGridFiles;
            end
        end

        function slices2frames(obj, prop, lims, label, area, name, slicePlane, runs)

            replot = false;

            if nargin < 5 || isempty(area)
                area = obj.blk.viewarea;
                aspect = obj.blk.aspect;
            else
                aspect = [(area(2)-area(1))...
                    (area(4)-area(3)) 1];
            end
            if nargin < 3 || isempty(lims)
                switch prop
                    case 'M'
                        lims = [0 1.4];
                    case 'tau_w'
                        lims = [0 800];
                    case 'vortZ'
                        lims = 1e5*[-0.7 0.7];
                    case 'T'
                        lims = [50 300];
                    case 'overlay'
                        lims = 5e4*[-1 1];
                end
            end
            if nargin < 4 || isempty(label)
                switch prop
                    case 'M'
                        label = 'M';
                    case 'tau_w'
                        label = '\tau_w';
                    case 'vortZ'
                        label = '\omega_z';
                    case 'T'
                        label = ('T (K)');
                end
            end

            if nargin < 7 || isempty(slicePlane)
                slicePlane = 'k';
            end

            if nargin < 8 || isempty(runs)
                runs = obj.run;
            end

            switch slicePlane
                case 'k'
                    var = ['k_' prop];
                case 'j'
                    var = ['j_' prop];
            end

            if isempty(runs)
                vname = [];
            elseif length(runs) == 1
                vname = ['/run' num2str(runs)];
            else
                vname = ['/run' num2str(runs(1)) '-' num2str(runs(end))];
            end
            if nargin < 6 || isempty(name)
                vname = [vname '_' var '.mp4"'];
            else
                vname = ['/run' num2str(runs(end)) '_' name '.mp4"'];
            end

            tempfolder = fullfile(obj.runpath,'temp');
            mkdir(tempfolder);

            [paths2read, slicenums2read, time2read, ishere] = obj.getSlicePaths(obj, runs);
            

            
            for i = 1:len(paths2read)

                imgfolder = fullfile(slices(i).casepath,'animation_images',var);
                if ~exist(imgfolder, 'dir')
                    mkdir(imgfolder);
                end
                fprintf('Plotting slice %d/%d\n',i,length(slices))
                fname = fullfile(imgfolder,sprintf('img_%03d.png',slices(i).nSlice));
                if ~exist(fname,'file') || replot

                
                    switch slicePlane
                        case 'k'
                            slice = kSlice(obj.blk,obj.gas,obj.bcs,paths2read{i},slicenums2read(i),time2read(i),obj.casetype,ishere);
    
                            repeats = 1;
                            if isfield(obj.blk, "n_pitchwise_repeats")
                                repeats = obj.blk.n_pitchwise_repeats;
                            end
                            h = figure('Visible','off');
                            ax = axes(h);
                          
                            q = slice.(prop);
                        
                            hold on
                            offset = 0;
                            if repeats > 2
                                offset = -obj.blk.pitch;
                            end
                            for ir = 1:repeats
                            for i=1:slice.NB
                                pcolor(ax, obj.blk.x{i}, obj.blk.y{i}+offset+(ir-1)*obj.blk.pitch, q{i});
                            end
                            end
                            shading('interp')
                            axis equal
                            pbaspect(aspect)
                            axis(area)
                            axis off
                            if ismember(string(prop),["vortZ","v","w"])
                                colormap(redblue)
                            end
                            cb = colorbar(ax);
                            if nargin > 4 && ~isempty(lims)
                                clim(lims);
                            end
                            if nargin > 5 && ~isempty(label)
                                cb.Label.String = label;
                            end
                            set(ax,'FontSize',12);
                            exportgraphics(h, fpath, 'Resolution', 600);
                        case 'j'
                            slice = jSlice(obj.blk,obj.gas,paths2read{i},slicenums2read(i),time2read(i),obj.casetype,ishere);
    
                    end
                end
                clear slice
            end

        end

        function frames2movie(obj,name,runs,framerate)
            if nargin < 3 || isempty(runs)
                runs = obj.runs;
            end
            if nargin < 4 || isempty(framerate)
                framerate = 5;
            end
            tempfolder = fullfile(obj.runpath,'temp');
            mkdir(tempfolder);
            if isempty(runs)
                folder = fullfile(obj.casepath,'animation_images',name);
                copyfile([folder '/*.png'], tempfolder);
            else
                for r = runs
                    folder = fullfile(obj.casepath,['run' num2str(r)], 'animation_images', name);
                    copyfile([folder '/*.png'], tempfolder);
                end
                vname = fullfile(obj.casepath, obj.casename, name);
            end
            if isempty(runs)
                vname = fullfile(obj.casepath, [obj.casename '_' name '.mp4']);
            elseif length(runs) == 1
                vname = fullfile(obj.casepath, ['run' num2str(runs(1)) '/' 'run' num2str(runs(1)) '_' name '.mp4']);
            else
                vname = fullfile(obj.casepath, ['run' num2str(runs(end)) '/' 'run' num2str(runs(1)) '-' ...
                   num2str(runs(end)) '_' name '.mp4']);
            end

            system(['ffmpeg -framerate ' num2str(framerate) ' -pattern_type glob -i "' tempfolder '/*.png" '...
                '-c:v libx264 -profile:v high -pix_fmt yuv420p -vf ' ...
                '"pad=ceil(iw/2)*2:ceil(ih/2)*2" "' ...
                 vname '"']);

            rmdir(tempfolder,"s");

        end

        function writeMovie(obj, slices, prop, varargin)
            

            defaultLims = [];
            defaultLabel = [];
            if ~isempty(obj.blk.viewarea)
                defaultViewArea = obj.blk.viewarea;
            else
                defaultViewArea = [];
            end
            if isfield(obj.blk, "n_pitchwise_repeats")
                defaultRepeats = obj.blk.n_pitchwise_repeats;
            else
                defaultRepeats = 1;
            end
            defaultRot = 0;

            p = inputParser;

            
            addRequired(p, 'slices', @slicesValidationFcn)
            addRequired(p, 'prop', @ischar);
            addParameter(p, 'lims', defaultLims);
            addParameter(p, 'label', defaultLabel);
            addParameter(p, 'viewarea', defaultViewArea);
            addParameter(p, 'nrepeats', defaultRepeats);
            addParameter(p, 'rot', defaultRot);
            addParameter(p, 'Interpreter', 'none');
            addParameter(p, 'ColorMap', 'parula')
            addParameter(p, 'name', []);
            addParameter(p, 'framerate', 10);
            addParameter(p, 'replot', false);

            function val = slicesValidationFcn(slices)
                
                if ischar(slices) || isstring(slices)
                    if ismember(string(slices), ["kSlices", "jSlices"])
                        val = true;
                    else
                        disp('Slices string must be "kSlices" or "jSlices"');
                        val = false;
                    end
                elseif ismember(string(class(slices)), ["kSlice", "jSlice"])
                    val = true;
                else
                    disp('Slices must be list of kSlices or jSlices');
                    val = false;
                end

            end

            parse(p, slices, prop, varargin{:});

            lims = p.Results.lims;
            label = p.Results.label;
            repeats = p.Results.nrepeats;
            framerate = p.Results.framerate;
            area = p.Results.viewarea;
            aspect = [(area(2)-area(1))...
                (area(4)-area(3)) 1];
            if isempty(lims)
                switch prop
                    case 'M'
                        lims = [0 1.4];
                    case 'tau_w'
                        lims = [0 800];
                    case 'vortZ'
                        lims = 1e5*[-0.7 0.7];
                    case 'T'
                        lims = [50 300];
                    case 'overlay'
                        lims = 5e4*[-1 1];
                    case 'w'
                        lims = [];
                end
            end
            if isempty(label)
                switch prop
                    case 'M'
                        label = 'M';
                    case 'tau_w'
                        label = '\tau_w';
                    case 'vortZ'
                        label = '\omega_z';
                    case 'T'
                        label = ('T (K)');
                end
            end

            readSlices = false;
            switch class(slices(1))
                case 'kSlice'
                    var = ['k_' prop];
                    sliceType = 'kSlice';
                case 'jSlice'
                    var = ['j_' prop];
                    sliceType = 'jSlice';
                case 'char'
                    readSlices = true;
                    [paths2read, slicenums2read, time2read, ishere] = obj.getSlicePaths(obj.run);
                    switch slices(1)
                        case 'k'
                            var = ['k_' prop];
                            sliceType = 'kSlice';
                        case 'j'
                            var = ['j_' prop];
                            sliceType = 'jSlice';
                    end
                    slices = slicenums2read;
            end

            if isempty(obj.run)
                vname = [];
            elseif length(obj.run) == 1
                vname = ['/run' num2str(obj.run)];
            else
                vname = ['/run' num2str(obj.run(1)) '-' num2str(obj.run(end))];
            end
            if isempty(p.Results.name)
                vname = [vname '_' var '.mp4"'];
            else
                vname = ['/run' num2str(obj.run(end)) '_' p.Results.name '.mp4"'];
            end

            tempfolder = fullfile(obj.runpath,'temp');
            mkdir(tempfolder);
            
            %%

            for i=1:length(slices)
                
                if readSlices
                    nSlice = slicenums2read(i);
                    imgfolder = fullfile(paths2read{i},'animation_images',var);
                else
                    nSlice = slices(i).nSlice;
                    imgfolder = fullfile(slices(i).casepath,'animation_images',var);
                end

                if ~exist(imgfolder, 'dir')
                    mkdir(imgfolder);
                end

                fprintf('Plotting slice %d/%d\n',i,length(slices))

                % Get path of frame and check if already plotted
                fname = fullfile(imgfolder,sprintf('img_%03d.png',nSlice));
                if ~exist(fname,'file') || p.Results.replot

                    % Get slice to contour
                    if readSlices
                        try
                            switch sliceType
                                case 'kSlice'
                                    slice = kSlice(obj.blk,obj.gas,obj.bcs,paths2read{i},slicenums2read(i),time2read(i),obj.casetype,ishere);
                                case 'jSlice'
                                    slice = jSlice(obj.blk,obj.gas,paths2read{i},slicenums2read(i),time2read(i),obj.casetype,ishere);
                            end
                        catch
                            continue
                        end
                    else 
                        slice = slices(i);
                    end

                    if strcmp(prop,'')%overlay')
                        slice2schlierenVortOverlay(slice, obj.blk, fullfile(imgfolder,sprintf('img_%03d.png',slice.nSlice)), lims, area, aspect);
                    else
                        switch class(slice)
                            case 'kSlice'
                                slice2kPlot(slice, obj.blk, prop, fullfile(imgfolder,sprintf('img_%03d.png',slice.nSlice)), lims, label, area, aspect);
                            case 'jSlice'
		                        slice2jPlot(slice, prop, fullfile(imgfolder, sprintf('img_%03d.png',slice.nSlice)), lims, label);
                        end
                    end

                    if readSlices
                        clear slice
                    end
                    
                end
                copyfile(fname, tempfolder);
            end
            
            system(['ffmpeg -framerate ' num2str(framerate) ' -pattern_type glob -i "' tempfolder '/*.png" '...
                '-c:v libx264 -profile:v high -pix_fmt yuv420p -vf ' ...
                '"pad=ceil(iw/2)*2:ceil(ih/2)*2" "' ...
                 obj.runpath vname]);

            rmdir(tempfolder,"s");

        end

        function [dim] = generateInflowTurb(obj, Tu, L, write)
            
            if nargin < 2
                write=true;
            end

            %L = 1.0;            % Lengthscale of largest eddies
            M = 10000;          % Number of Fourier modes
            %Tu = 0.005;         % Turbulence intensity

            [xb,yb,zb,dx] = writeTurbGrid(obj.blk.inlet_blocks{1}, obj.solver.span, obj.casepath);
            [lMin, lMax, minL, ~] = gridSpac(xb,yb,zb);
            [ni,nj,nk] = size(xb);
            dim = [ni,nj,nk];

            rgas = obj.gas.cp*(obj.gas.gam-1)/obj.gas.gam;
            Tin = obj.bcs.Toin - obj.bcs.vin^2/(2*obj.gas.cp);
            ronow = (obj.bcs.Poin/(rgas*Tin)) * (Tin/obj.bcs.Toin)^(obj.gas.gam/(obj.gas.gam-1));
            nunow = 1/ronow*obj.gas.mu_ref*((obj.gas.mu_tref+obj.gas.mu_cref)/(Tin+obj.gas.mu_cref))*(Tin/obj.gas.mu_tref)^1.5;

            if write
                dir = pwd;
                cd(fullfile(obj.casepath,'turb-gen'))
                f = fopen('turbin_params_grid.txt','w');
                fprintf(f,'%s\n','turb-grid.dat');
                fprintf(f,'%d, %d, %d\n',[ni nj nk]);
                fprintf(f,'%0.6e\n', lMax);
                fprintf(f,'%e\n', minL);
                fprintf(f,'%d\n',M);
                fprintf(f,'%e\n', nunow);
                fprintf(f,'%f, %f\n', [Tu obj.bcs.vin]);
                fprintf(f,'%f, %e\n', [-1 L]);
                fprintf(f,'%s\n', 'inflow_turb_gridTurb.dat');
                fclose(f)
                system('turbInGrid > out &')
                cd(dir)
            end

            vIn = obj.bcs.vin;
            eps =  (Tu*vIn)^3/L*sqrt(3/2)^3;
            nu = nunow;
            tol = 0.00001;
            wrap = @(ke,keta,dltk,b,tol) estInt(ke,keta,dltk,b,tol); 

            %rms of velocity fluctuations
            u = Tu*vIn;
            %ch'ic length scale-based Re
            %ReL = u*Lc/nui;
            ReL = u*L/(nu)
            %Wave length related to location of energy maximum (at sqrt(12/5)ke)
            %relationship to int. length scale fixed for isotropic turbulence
            keSt = 0.746834/(L*0.09);
            %Kolmogorov wave length (= eps^1/4/nu^3/4 = u^(3/4)/(nu^(3/4)*Lint^(1/4))
            %keta = ReL.^(3/4)/(Lint^(1/4)*Lc^(3/4));
            keta = eps^(1/4)/nu^(3/4);
        
            %Grid spacing
            %dg = maxL/(N-1)
            %minimum resolved wave number (i.e. longest waves - want to take minL
            %here rather than maxL, so that we won't end up with waves too large
            %for one direction.
            kmi = 2*pi/minL;
            %maximum resolved wave number
            kmx = pi/lMin;  
        
        
            mn = linspace(1,M,M);  
            km = zeros(M,1);
        
            ce = 2*nu;
            cu = 2/3;
            b = ce*u^2/(eps*cu);
        
            %set wavenumber as centre of bin
            for m = 1:M    
                km(m) = kmi+(kmx-kmi)/M*(mn(m)-0.5);
            end
        
            %width of wave-number bins
            dltk = (max(km)-min(km))/(M-1);
        
            try
                %Find ke..
                fun = @(x) wrap(x,keta,dltk,b,tol);    % function of x alone
            catch
                
            end
            ke = fzero(fun,[1, keSt]);
            
            disp(['Calculated k_e as ', num2str(ke)])
            [dsint, kint, alpha] = calcAlpha(ke,dltk,keta,nu,u,tol);
        
            
            disp(['Calculated alpha as ', num2str(alpha)])
        
            Ek = zeros(M,1);
        
            for m = 1:M
                bfk = km(m)/ke;
                expc = exp(-2*(km(m)/keta)^2);
                Ek(m) = alpha*u^2/ke*bfk^4/(1+bfk^2)^(17/6)*expc;
            end
        
            E = Ek/u^2*keta;
            k = km/keta;
            
            disp(1/u^2*keta)
            disp(keta)
            
            figure()
            
            loglog(k,E,'LineWidth',1.5);
            [~,minDex] = min(abs(k-1.5));
            xLimits = [min(k) 1.5];                   %# Limits for the x axis
            yLimits = [E(minDex) max(E)*1.1 ];                      %# Limits for the y axis
            %logScale = diff(yLimits)/diff(xLimits);  %# Scale between the x and y ranges
            %powerScale = diff(log10(yLimits))/...    %# Scale between the x and y powers
            %             diff(log10(xLimits));
            %set(gca,'Xlim',xLimits,'YLim',yLimits,...              %# Set the limits and the
            %        'DataAspectRatio',[1 logScale/powerScale 1]);  %#   data aspect ratio
            set(gca,'Xlim',xLimits,'YLim',yLimits);    
            grid on
            hold on
            i1 = find(E>=0.8*max(E),1,'last');
            i2 = find(E<=0.05*max(E),1,'first');
            Bp = polyfit(log10(k(i1:i2)), log10(E(i1:i2)), 1);
            Yp = polyval(Bp,log10(k(i1:i2)));
            Yp=10.^Yp;
            loglog(k(i1:i2),Yp,'-r');
            
            %find approximate gradient in 
            m = (log10(Yp(end))-log10(Yp(1)))/(log10(k(i2))-log10(k(i1)))
            
            keta
            xlabel('$\kappa/\kappa_\eta$','interpreter','latex','FontSize',15)
            ylabel('$E(\kappa)\kappa_\eta/u''^2$','interpreter','latex','FontSize',15)
            
            fprintf('Set ilength to %d\n',ni);
            fprintf('Set lturb to: %9.7f\n', dx);

            obj.bcs.ilength = ni;
            obj.bcs.lturb = dx;


        end

        function [u,v,w] = conditionInflowTurb(obj,dim)
            gridPath = py.str(fullfile(obj.casepath,'turb-gen','turb-grid.dat'));
            namePath = py.str(fullfile(obj.casepath,'turb-gen','inflow_turb_gridTurb.dat'));
            DIM = py.tuple(int32(dim));
            np = py.importlib.import_module('numpy');
            turb_postproc = py.importlib.import_module('turb_postproc');
            U = turb_postproc.conditionAxSlc(namePath,gridPath,DIM,true,false,false);
            u = U(1); v = U(2); w = U(3);
            namePath = py.str(fullfile(obj.casepath,'turb-gen','inflow_turb.dat'));

            u = double(np.ascontiguousarray(u{1}));
            v = double(np.ascontiguousarray(v{1}));
            w = double(np.ascontiguousarray(w{1}));

            turb_postproc.write3DNSorder(namePath,u,v,w)

        end

        function [u,v,w] = readTurb3DNSOrder(obj,dim)
            gridPath = py.str(fullfile(obj.casepath,'turb-gen','turb-grid.dat'));
            namePath = py.str(fullfile(obj.casepath,'turb-gen','inflow_turb.dat'));
            DIM = py.tuple(int32(dim));
            np = py.importlib.import_module('numpy');
            turb_postproc = py.importlib.import_module('turb_postproc');


            flo = turb_postproc.read3DNSorder(namePath,DIM);
            u = double(np.ascontiguousarray(flo{1}));
            v = double(np.ascontiguousarray(flo{2}));
            w = double(np.ascontiguousarray(flo{3}));

        end



        function write_HYDRA_mesh(obj)
            path = fullfile(obj.casepath,[obj.casename '_HYDRA.xyz']);
            nk = 6;
            span = (nk-1)*obj.solver.span/(obj.solver.nk-1);
            write_plot3d_extruded(obj.blk, path, nk, span);
        end

        function blkNodes = writeFluentMesh2d(obj, path)
            if nargin < 2
                path = fullfile(obj.casepath,[obj.casename '_Fluent_2d.msh']);
                %npath = fullfile(obj.casepath,[obj.casename '_Fluent_nodes.mat']);
            end
            bnd = obj.getBoundaries;
            blkNodes = writeFluentMesh(path, obj.blk, obj.blk.next_block, obj.blk.next_patch, bnd, true);
        end

        function blkNodes = writeFluentMeshExtruded(obj, spannow, nknow, path)
            bnd = obj.getBoundaries;
            if nargin < 4
                path = fullfile(obj.casepath,[obj.casename '_extruded.msh']);
            end
            blkNodes = writeFluentMeshExtruded(path, obj.blk, bnd, spannow, nknow, true);
        end


%         function blkNodes = write_fluent_mesh_2d(obj, path)
%             blkNodes = writeAerofoilFluentMesh(path, obj.blk, obj.blk.next_block, obj.blk.next_patch, true);
%         end

        function blkNodes = write_Fluent_mesh(obj, fname)
            if nargin < 2
                path = fullfile(obj.casepath,[obj.casename '_Fluent.msh']);
                %npath = fullfile(obj.casepath,[obj.casename '_Fluent_nodes.mat']);
            else
                path = fullfile(obj.casepath,[fname '.msh']);
            end

            nk = 6;
            span = (nk-1)*obj.solver.span/(obj.solver.nk-1);
            blkNodes = writeFluentMeshExtruded(path, obj.blk, obj.blk.next_block, obj.next_patch, nk, span, true);
        end

        function write_2d_plot3d_mesh(obj)
            path = fullfile(obj.casepath,[obj.casename '_2d.xyz']);
            write_plot3d_2d(obj.blk, path);
        end

        function writeFluentCas2D(obj, path)
            if nargin < 2
                path = fullfile(obj.casepath, [obj.casename '_2d.cas']);
            end
            writeFluentCas(path, obj.blk, obj.getBoundaries);
        end

        function writeFluentCas3D(obj, path)
            if nargin < 2
                path = fullfile(obj.casepath, [obj.casename '.cas']);
            end
            writeFluentCas3D(path, obj.blk, obj.getBoundaries);
        end

        function slice = vol2jSlice(obj)

            slice = jSlice(obj.blk, obj.gas, obj.bcs);


            x = [];
            xs = [];
            y = [];
            y0 = [];

            for iblk=1:length(obj.blk.oblocks)

                ib = obj.blk.oblocks(iblk);
                
                xtmp = obj.blk.x{ib}(:,end-1);
                xstmp = obj.blk.x{ib}(:,end);
                ytmp = obj.blk.y{ib}(:,end-1);
                dxtmp = obj.blk.x{ib}(:,end-1) - obj.blk.x{ib}(:,end);
                dytmp = obj.blk.y{ib}(:,end-1) - obj.blk.y{ib}(:,end);
                y0tmp = sqrt(dxtmp.^2+dytmp.^2);

                ro = obj.instFlow.ro{ib}(:,end-1,:);
                u = obj.instFlow.u{ib}(:,end-1,:);
                v = obj.instFlow.v{ib}(:,end-1,:);
                w = obj.instFlow.w{ib}(:,end-1,:);
                Et = obj.instFlow.Et{ib}(:,end-1,:);
                
                if obj.blk.oblocks_flip(iblk) == 1
                    ro = flip(ro);
                    u = flip(u);
                    v = flip(v);
                    w = flip(w);
                    Et = flip(Et);
                    xtmp = flip(xtmp);
                    ytmp = flip(ytmp);
                    y0tmp = flip(y0tmp);
                    xstmp = flip(xstmp);
                end
                
                slice.ro = [slice.ro; ro(1:end-1,:)];
                slice.u = [slice.u; u(1:end-1,:)];
                slice.v = [slice.v; v(1:end-1,:)];
                slice.w = [slice.w; w(1:end-1,:)];
                slice.Et = [slice.Et; Et(1:end-1,:)];
                x = [x; xtmp(1:end-1)];
                y = [y; ytmp(1:end-1)];
                y0 = [y0; y0tmp(1:end-1)];
                xs = [xs; xstmp(1:end-1)];

            end
            x = x';
            y = y';
            y0 = y0';
            [~, iLE] = min(xs);
            [~, iTE] = max(xs);
            slice.ro = slice.ro(iLE:iTE,:);
            slice.u = slice.u(iLE:iTE,:);
            slice.v = slice.v(iLE:iTE,:);
            slice.w = slice.w(iLE:iTE,:);
            slice.Et = slice.Et(iLE:iTE,:);

            x = x(iLE:iTE);
            y = y(iLE:iTE);
            y0 = y0(iLE:iTE);
            slice.ssurf(1) = 0;
            for i=2:length(x)
                dx = x(i)-x(i-1);
                dy = y(i)-y(i-1);
                ds = sqrt(dx^2 + dy^2);
                slice.ssurf(i) = slice.ssurf(i-1)+ds;
            end
            z = linspace(0,obj.blk.span,obj.blk.nk);
            [slice.Z, slice.X] = meshgrid(z,slice.ssurf);
            slice.Y0 = repmat(y0',1,obj.blk.nk);
            
        end

        function mean = inst2ave(obj, sliceNums)
            mean = meanSlice(obj.blk, obj.gas, obj.bcs, []);

            for ib=1:obj.NB
                ni = obj.blk.blockdims(ib,1);
                nj = obj.blk.blockdims(ib,2);

                ro{ib} = zeros(ni,nj);
                ru{ib} = zeros(ni,nj);
                rv{ib} = zeros(ni,nj);
                rw{ib} = zeros(ni,nj);
                Et{ib} = zeros(ni,nj);

                u{ib} = zeros(ni, nj);
                v{ib} = zeros(ni, nj);
                w{ib} = zeros(ni, nj);

                uu{ib} = zeros(ni,nj);
                vv{ib} = zeros(ni,nj);
                ww{ib} = zeros(ni,nj);

                ro2{ib} = zeros(ni,nj);
                rou2{ib} = zeros(ni,nj);
                rov2{ib} = zeros(ni,nj);
                row2{ib} = zeros(ni,nj);

                rouv{ib} = zeros(ni,nj);
                rouw{ib} = zeros(ni,nj);
                rovw{ib} = zeros(ni,nj);

                pbar{ib} = zeros(ni,nj);
                Tbar{ib} = zeros(ni,nj);

                roUddUdd{ib} = zeros(ni, nj);
                roVddVdd{ib} = zeros(ni, nj);
                roWddWdd{ib} = zeros(ni, nj);

                roUddUdd{ib} = zeros(ni, nj);
                roUddWdd{ib} = zeros(ni, nj);
                roVddWdd{ib} = zeros(ni, nj);

                Pr{ib} = zeros(ni, nj);
                k{ib} = zeros(ni, nj);

            end
            
            [paths2read, slicenums2read, time2read, ishere] = obj.getSlicePaths(obj.run);
            if nargin > 1 && ~isempty(sliceNums)
                [~, inds] = ismember(slicenums2read, sliceNums);
                inds = inds(inds>0);
                nSlices = length(inds);
            else
                nSlices = length(slicenums2read);
                inds = 1:nSlices;
            end

            for i = 1:nSlices

                ind = inds(i);

                fprintf('Reading slice %d/%d\n', [i, nSlices])
                slice = kSlice(obj.blk, obj.gas, obj.bcs, paths2read{ind}, slicenums2read(ind), time2read(ind), obj.casetype, ishere);
                
                for ib = 1:obj.NB
                    ro{ib} = ro{ib} + slice.ro{ib} / nSlices;
                    ru{ib} = ru{ib} + slice.ro{ib} .* slice.u{ib} / nSlices;
                    rv{ib} = rv{ib} + slice.ro{ib} .* slice.v{ib} / nSlices;
                    rw{ib} = rw{ib} + slice.ro{ib} .* slice.w{ib} / nSlices;
                    Et{ib} = Et{ib} + slice.Et{ib} / nSlices;

                    u{ib} = u{ib} + slice.u{ib} / nSlices;
                    v{ib} = v{ib} + slice.v{ib} / nSlices;
                    w{ib} = w{ib} + slice.w{ib} / nSlices;

                    uu{ib} = uu{ib} + slice.u{ib}.^2 / nSlices;
                    vv{ib} = vv{ib} + slice.v{ib}.^2 / nSlices;
                    ww{ib} = ww{ib} + slice.w{ib}.^2 / nSlices;

                    ro2{ib} = ro2{ib} + slice.ro{ib}.^2 / nSlices;
                    rou2{ib} = rou2{ib} + slice.ro{ib} .* slice.u{ib}.^2 / nSlices;
                    rov2{ib} = rov2{ib} + slice.ro{ib} .* slice.v{ib}.^2 / nSlices;
                    row2{ib} = row2{ib} + slice.ro{ib} .* slice.w{ib}.^2 / nSlices;

                    rouv{ib} = rouv{ib} + slice.ro{ib} .* slice.u{ib} .* slice.v{ib} / nSlices;
                    rouw{ib} = rouw{ib} + slice.ro{ib} .* slice.u{ib} .* slice.w{ib} / nSlices;
                    rovw{ib} = rovw{ib} + slice.ro{ib} .* slice.v{ib} .* slice.w{ib} / nSlices;

                    
                end

                clear slice

            end
            
            for ib = 1:obj.NB
                % u{ib} = u{ib} ./ ro{ib};
                % v{ib} = v{ib} ./ ro{ib};
                % w{ib} = w{ib} ./ ro{ib};
                pbar{ib} = (Et{ib} - 0.5*(rou2{ib} + rov2{ib} + row2{ib}))*(obj.gas.gam-1);
                Tbar{ib} = (pbar{ib}.*obj.gas.gam)./(obj.gas.cp*(obj.gas.gam-1)*ro{ib});


                [DUDX,DUDY] = gradHO(obj.blk.x{ib},obj.blk.y{ib},u{ib});
                [DVDX,DVDY] = gradHO(obj.blk.x{ib},obj.blk.y{ib},v{ib});

                % UdUd = rou2{ib}./ro{ib} - u{ib}.*u{ib};
                % VdVd = rov2{ib}./ro{ib} - v{ib}.*v{ib};
                % WdWd = row2{ib}./ro{ib} - w{ib}.*w{ib};

                % UdVd = rouv{ib}./ro{ib} - u{ib}.*v{ib};
                % UdWd = rouw{ib}./ro{ib} - u{ib}.*w{ib};
                % VdWd = rovw{ib}./ro{ib} - v{ib}.*w{ib};

                roUddUdd{ib} = rou2{ib} - ro{ib}.*u{ib}.*u{ib};
                roVddVdd{ib} = rov2{ib} - ro{ib}.*v{ib}.*v{ib};
                roWddWdd{ib} = row2{ib} - ro{ib}.*w{ib}.*w{ib};

                roUddVdd{ib} = rouv{ib} - ro{ib}.*u{ib}.*v{ib};
                roUddWdd{ib} = rouw{ib} - ro{ib}.*u{ib}.*w{ib};
                roVddWdd{ib} = rovw{ib} - ro{ib}.*v{ib}.*w{ib};

                % Pr{ib} = -ro{ib}.*(UdUd.*DUDX + UdVd.*(DUDY+DVDX) + VdVd.*DVDY);
                % k{ib} = 0.5*(UdUd + VdVd + WdWd);

                UdUd{ib} = uu{ib} - u{ib}.*u{ib};
                VdVd{ib} = vv{ib} - v{ib}.*v{ib};
                WdWd{ib} = ww{ib} - w{ib}.*w{ib};

                RdUd{ib} = ru{ib} - ro{ib}.*u{ib};
                RdVd{ib} = rv{ib} - ro{ib}.*v{ib};
                RdWd{ib} = rw{ib} - ro{ib}.*w{ib};

                RUdUd{ib} = rou2{ib} - ro{ib}.*u{ib}.*u{ib} - 2*u{ib}.*RdUd{ib};
                RVdVd{ib} = rov2{ib} - ro{ib}.*v{ib}.*v{ib} - 2*v{ib}.*RdVd{ib};
                
                RUdVd{ib} = rouv{ib} - ro{ib}.*u{ib}.*v{ib} - u{ib}.*RdVd{ib} - v{ib}.*RdUd{ib};

                Pr{ib} = -1 * (RUdUd{ib}.*DUDX + RUdVd{ib}.*(DUDY+DVDX) + RVdVd{ib}.*DVDY);
                k{ib} = 0.5 * (UdUd{ib} + VdVd{ib} + WdWd{ib});

            end


            mean.ro = ro;
            mean.u = u;
            mean.v = v;
            mean.w = w;
            mean.Et = Et;
    
            mean.pbar = pbar;
            mean.Tbar = Tbar;
    
            mean.roUddUdd = roUddUdd;
            mean.roVddVdd = roVddVdd;
            mean.roWddWdd = roWddWdd;
    
            mean.roUddVdd = roUddVdd;
            mean.roUddWdd = roUddWdd;
            mean.roVddWdd = roVddWdd;
    
            mean.Pr = Pr;
            mean.k = k;

            obj.meanFlow = mean;

        end



        function update_Min(obj, M)
            obj.bcs.vin = Vel_M(M, obj.bcs.Toin, obj.gas.cp, obj.gas.gam);
            obj.bcs.pexit = obj.bcs.Poin*p_p0(M, obj.gas.gam);
            obj.writeInputFiles;
        end

    end
end
