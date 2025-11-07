classdef kStatsSlice < kCut
    % KSLICE Contains a 2D slice of the flow at k-boundary
    %   Detailed explanation goes here

    properties
        time;
        nSlice;
        casepath;
    end

    properties (Dependent = true, Hidden = true)
%         T;
%         p;
%         M;
%         s;
%         vel;
    end

    methods
        function obj = kStatsSlice(blk, gas, bcs, casedir, nSlice, time, casetype, ishere)
            obj@kCut(blk, gas, bcs);
            disp('Constructing kStatsSlice')

            if nargin < 8
                ishere = false;
            end
        
            if nargin > 3
    
                if nargin > 2
                    if ~ischar(nSlice)
                        obj.nSlice = nSlice;
                        obj.time = time;
                    end

                    obj.casepath = casedir;
        
                    for nb = 1:obj.NB

                        ni = blk.blockdims(nb,1);
                        nj = blk.blockdims(nb,2);

                        ro = zeros(ni,nj);
                        ru = zeros(ni,nj);
                        rv = zeros(ni,nj);
                        rw = zeros(ni,nj);
                        Et = zeros(ni,nj);

                        switch casetype
                            case 'cpu'
                                % if ishere
                                %     flopath = fullfile(casedir,  ['kcu2_' num2str(nb) '_' num2str(nSlice)]);
                                %     flofile = fopen(flopath,'r');
                                %     nodfile = fopen(fullfile(casedir, ['knd2_' num2str(nb) '_' num2str(nSlice)]),'r');
                                % else
                                %     flopath = fullfile(casedir, 'k_cuts',  ['kcu2_' num2str(nb) '_' num2str(nSlice)]);
                                %     flofile = fopen(flopath,'r');
                                %     nodfile = fopen(fullfile(casedir, 'k_cuts', ['knd2_' num2str(nb) '_' num2str(nSlice)]),'r');
                                % end
                                % A = fread(flofile,inf,'float64');
                                % A = reshape(A,5,length(A)/5);
                                % 
                                % B = fread(nodfile,inf,'uint32');
                                % B = reshape(B,3,length(B)/3);
                                % 
                                % fclose(flofile);
                                % fclose(nodfile);
                                % 
                                % for n=1:size(A,2)
                                %     i = B(1,n);
                                %     j = B(2,n);
                                %     k = B(3,n);
                                %     ro(i,j) = A(1,n);
                                %     ru(i,j) = A(2,n);
                                %     rv(i,j) = A(3,n);
                                %     rw(i,j) = A(4,n);
                                %     Et(i,j) = A(5,n);
                                % end

                            case 'gpu'
                                if ischar(nSlice)
                                    fid = fopen(fullfile(casedir, [nSlice '_' num2str(nb)]));
                                else
                                    if ishere
                                        fid = fopen(fullfile(casedir, ['kcut_' num2str(nb) '_' num2str(nSlice)]));
                                    else
                                        fid = fopen(fullfile(casedir, 'k_cuts', ['kbar_' num2str(nb) '_' num2str(nSlice)]));
                                    end
                                end
                                nstats_prim = 11;
                                A = fread(fid, ni*nj*nstats_prim, 'float64');
                                A = reshape(A, nstats_prim, length(A)/nstats_prim)';
                                
                                nstats_entr = 3;
                                B = fread(fid, ni*nj*nstats_entr, 'float64');
                                B = reshape(B, nstats_entr, length(B)/nstats_entr)';

                                ro = reshape(A(:,1),ni,nj);
                                ru = reshape(A(:,2),ni,nj);
                                rv = reshape(A(:,3),ni,nj);
                                rw = reshape(A(:,4),ni,nj);
                                Et = reshape(A(:,5),ni,nj);

                                % ruu = reshape(A(:,6),ni,nj);
                                % rvv = reshape(A(:,7),ni,nj);
                                % rww = reshape(A(:,8),ni,nj);
                                % ruv = reshape(A(:,9),ni,nj);
                                % ruw = reshape(A(:,10),ni,nj);
                                % rvw = reshape(A(:,11),ni,nj);
                                % rus = reshape(A(:,12),ni,nj);
                                % rvs = reshape(A(:,13),ni,nj);
                                % lnT = reshape(A(:,614),ni,nj);
                                
                        end
        
                        obj.ro{nb} = ro;
                        obj.u{nb} = ru./ro;
                        obj.v{nb} = rv./ro;
                        obj.w{nb} = rw./ro;
                        obj.Et{nb} = Et;
                    end
                end
            end
        end



%         function value = get.vel(obj)
%             value = cell(1,obj.NB);
%             for nb =1:obj.NB
%                 value{nb} = sqrt(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2);
%             end
%         end
% 
%         function value = get.M(obj)
%             value = cell(1,obj.NB);
%             for nb =1:obj.NB
%                 value{nb} = obj.vel{nb}./sqrt(obj.gas.gam*obj.gas.rgas*obj.T{nb});
%             end
%         end
% 
%         function value = get.s(obj)
%             value = cell(1,obj.NB);
%             for nb =1:obj.NB
%                 value{nb} = obj.gas.cp*log(obj.T{nb}/300) - obj.gas.rgas*log(obj.p{nb}/1e5);
%             end
%         end

%         function value = get.vortZ(obj)
%             value = cell(1,obj.NB);
%             for nb =1:obj.NB
%                 [~,DUDY] = gradHO(obj.blk.x{nb},obj.blk.y{nb},obj.u{nb});
%                 [DVDX,~] = gradHO(obj.blk.x{nb},obj.blk.y{nb},obj.v{nb});
%                 value{nb} = DVDX-DUDY;
%             end
%         end

        function getSize(obj)
            props = properties(obj); 
            totSize = 0; 
           
            for ii=1:length(props) 
                currentProperty = getfield(obj, char(props(ii))); 
                temp = whos('currentProperty'); 
                totSize = totSize + temp.bytes; 
            end
          
            fprintf(1, '%d MB\n', totSize/1e6);
        end

        function value = inst2ave(obj)
            value = aveSlice(obj.blk, obj.gas);
            value.ro = obj.ro;
            value.u = obj.u;
            value.v = obj.v;
            value.w = obj.w;
            value.Et = obj.Et;
        end
    end
end