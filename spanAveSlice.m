classdef spanAveSlice < aveSlice
    %SPANAVESLICE Class to store a spanwise averaged instantaneous flow

    properties
        nMean;
        ros_store = [];
    end

    methods
        function obj = spanAveSlice(blk,gas,bcs,casedir,casetype,nMean)
            obj@aveSlice(blk, gas, bcs);
            disp('Constructing spanAveSlice')
            
            if nargin > 3
                
                obj.nMean = nMean;

                for nb = 1:obj.NB
                    ni = blk.blockdims(nb,1);
                    nj = blk.blockdims(nb,2);
    
                    ro = zeros(ni,nj);
                    ru = zeros(ni,nj);
                    rv = zeros(ni,nj);
                    rw = zeros(ni,nj);
                    Et = zeros(ni,nj);
                    ros = zeros(ni,nj);

                    switch casetype
                        case 'cpu'
                            
                        case 'gpu'
                            nstats = 6;
                            flowpath = fullfile(casedir, ['flow2_' num2str(nb) '_' num2str(obj.nMean)]);
                            f = dir(flowpath);
                            flofile = fopen(flowpath);
                            nstats = f.bytes/(ni*nj*8);
                            A = fread(flofile,ni*nj*nstats,'float64');
                            A = reshape(A,nstats,[])';
                            fclose(flofile);

                            ro = reshape(A(:,1),ni,nj);
                            ru = reshape(A(:,2),ni,nj);
                            rv = reshape(A(:,3),ni,nj);
                            rw = reshape(A(:,4),ni,nj);
                            Et = reshape(A(:,5),ni,nj);
                            if nstats > 5
                                ros = reshape(A(:,6),ni,nj);
                            end
                    end


                    u = ru./ro;
                    v = rv./ro;
                    w = rw./ro;

                    obj.ro{nb} = ro;
                    obj.u{nb} = u;
                    obj.v{nb} = v;
                    obj.w{nb} = w;
                    obj.Et{nb} = Et;

                    if nstats > 5
                        obj.ros_store{nb} = ros;
                    end

                end
            end
        end

        function value = get_ros(obj)
            if ~isempty(obj.ros_store)
                value = obj.ros_store;
            else
                value = cell(1,obj.NB);
                Tnow = obj.T;
                pnow = obj.p;
                for nb = 1:obj.NB
                    value{nb} = obj.ro{nb}.*(obj.gas.cp*log(Tnow{nb}/obj.bcs.Toin) - ...
                        obj.gas.cp*(1-1/obj.gas.gam)*log(pnow{nb}/obj.bcs.Poin));
                end
            end
        end
    end
end