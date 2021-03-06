classdef meanSlice < aveSlice
    % MEANSLICE Contains the 2D (spanwise averaged) mean flow

        
    properties
        time;
        nMean;
        meanTime;
        Pr;             % Turbulence production
        diss;
        rev_gen_x;
        rev_gen_y;
        irrev_gen;
        diss_T;
        p;
        T;
        span;
        nk;
        rous;
        rovs;
        k;
        advK;
    end

    properties (Dependent = true)

    end

    methods
        function obj = meanSlice(casedir, blk, gas)
            obj@aveSlice(blk, gas);
            disp('Constructing meanSlice')
            if exist(fullfile(casedir,'mean_flo','nstats.txt'),'file')
                fid = fopen(fullfile(casedir,'mean_flo','nstats.txt'));
                nstats = str2double(fgetl(fid));
                statstype = str2double(fgetl(fid));
                fclose(fid);
            else
                nstats = 17;
                statstype = 1;
            end
            nstats

            if nargin > 0

                fullfile(casedir, 'mean_flo', 'mean_time.txt');
                fid = fopen(fullfile(casedir, 'mean_flo', 'mean_time.txt'));
                while ~feof(fid) % Use lastest mean files
                    temp=fgetl(fid);
                end
%                 temp = fgetl(fid);
                fclose(fid);
                temp = str2num(temp);
                obj.nMean = temp(1);
                obj.meanTime = temp(3);
                for nb = 1:obj.NB

                    flopath = fullfile(casedir, 'mean_flo',  ['mean2_' num2str(nb) '_' num2str(obj.nMean)]);
                    flofile = fopen(flopath,'r');
                    fullfile(casedir, 'mean_flo', ['mnod2_' num2str(nb) '_' num2str(obj.nMean)]);
                    nodfile = fopen(fullfile(casedir, 'mean_flo', ['mnod2_' num2str(nb) '_' num2str(obj.nMean)]),'r');
                    A = fread(flofile,inf,'float64');
                    A = reshape(A,nstats,length(A)/nstats);
                    
                    B = fread(nodfile,inf,'uint32');
                    B = reshape(B,3,length(B)/3);
            
                    fclose(flofile);
                    fclose(nodfile);
    
                    rodt = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    rudt = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    rvdt = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    rwdt = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    Etdt = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));

                    ro2 = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    rou2 = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    rov2 = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    row2 = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));

                    rouv = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    rouw = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    rovw = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));

                    p2 = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    p = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    T = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    rous = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    rovs = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    rows = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    diss = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    qx_T = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    qy_T = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    qz_T = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    irrev_gen = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));
                    diss_T = zeros(blk.blockdims(nb,1),blk.blockdims(nb,2));

                    sz = size(ro2);
                    icount(1:prod(sz)) = 0;
                    for n=1:size(A,2)

                        i = B(1,n);
                        j = B(2,n);

                        if icount(sub2ind(sz,i,j)) == 0

                            icount(sub2ind(sz,i,j)) = 1;

                            rodt(i,j) = A(1,n);
                            rudt(i,j) = A(2,n);
                            rvdt(i,j) = A(3,n);
                            rwdt(i,j) = A(4,n);
                            Etdt(i,j) = A(5,n);
    
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
                                irrev_gen(i,j) = A(21,n)/obj.meanTime;
                            end

                            if nstats > 21
                                diss_T(i,j) = A(22,n)/obj.meanTime;
                            end
                        end
                    end
    
                    obj.ro{nb} = rodt/obj.meanTime;
                    obj.u{nb} = rudt./(rodt);
                    obj.v{nb} = rvdt./(rodt);
                    obj.w{nb} = rwdt./(rodt);
                    obj.Et{nb} = Etdt/obj.meanTime;
                    obj.p{nb} = (Etdt/obj.meanTime - 0.5*(rou2 + rov2 + row2))*(obj.gas.gam-1);
                    obj.T{nb} = (obj.p{nb}.*obj.gas.gam)./(obj.gas.cp*(obj.gas.gam-1)*obj.ro{nb});
                    obj.diss{nb} = diss;
                    obj.rous{nb} = rous;
                    obj.rovs{nb} = rovs;

                    [DUDX,DUDY] = gradHO(blk.x{nb},blk.y{nb},obj.u{nb});
                    [DVDX,DVDY] = gradHO(blk.x{nb},blk.y{nb},obj.v{nb});

                    UdUd = rou2./obj.ro{nb} - obj.u{nb}.*obj.u{nb};
                    UdVd = rouv./obj.ro{nb} - obj.u{nb}.*obj.v{nb};
                    VdVd = rov2./obj.ro{nb} - obj.v{nb}.*obj.v{nb};
                    WdWd = row2./obj.ro{nb} - obj.w{nb}.*obj.w{nb};

                    obj.Pr{nb} = -(UdUd.*DUDX + UdVd.*(DUDY+DVDX) + VdVd.*DVDY);
                    obj.diss{nb} = diss;
                    obj.k{nb} = 0.5*(UdUd + VdVd + WdWd);

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
                
            end
            obj.span = blk.span;
            obj.nk = blk.nk{1};
            obj.getBCs(blk.inlet_blocks{1});
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

        
        

%         function value = get.dsdy(obj)
%             s = obj.oGridProp('s');
%             value = (s(:,2:end)-s(:,1:end-1))./(obj.yBL(:,2:end)-obj.yBL(:,1:end-1));
%         end
%         
% 
%         function value = get.Msurf(obj)
%             disp('Calculating surface M')
%             psurf = [];
%             pnow = obj.p;
%             size(pnow{4})
%             for i=1:length(obj.oblocks)
%                 clear temp
%                 temp = pnow{obj.oblocks(i)}(:,end);
%                 size(temp)
%                 if obj.oblocks_flip(i) == 1
%                     temp = flip(temp);
%                 end
%                 psurf = [psurf temp'];
%             end
%             value = sqrt((2/(obj.gas.gam - 1)) * ( (psurf/obj.p0in).^(-(obj.gas.gam-1)/obj.gas.gam) - 1));
%         end
% 
%         function value = get.BLedgeInd(obj)
%             temp = obj.dsdy;
%             for i=1:size(temp,1)
%                 j=3;
%                 temp(i,j);
%                 while temp(i,j) < obj.dsdyThresh
%                     j = j+1;
%                 end
%                 value(i) = j+1;
%             end
%         end
% 
%         function value = get.delta99(obj)
%             inds = obj.BLedgeInd;
%             for i=1:size(obj.yBL,1)
%                 value(i) = obj.yBL(i,inds(i));
%             end
%         end
% 
%         function value = get.U(obj)
%             disp('Calculating U')
%             unow = obj.oGridProp('u');
%             vnow = obj.oGridProp('v');
%             nnow = obj.n;
%             value = zeros(size(obj.yBL));
%             for i=1:size(obj.yBL,1)
%                 for j=1:size(obj.yBL,2)
%                     velnow = [unow(i,j); vnow(i,j)];
%                     value(i,j) = norm(velnow - nnow(:,i)*dot(nnow(:,i),velnow));
%                 end
%             end
%         end
% 
%         function value = get.delStar(obj)
%             inds = obj.BLedgeInd;
%             ronow = obj.oGridProp('ro');
%             Unow = obj.U;
%             value = zeros(1,length(inds));
%             for i=1:size(obj.yBL,1)
%                 integrand = 1 - ronow(i,1:inds(i)).*Unow(i,1:inds(i))/(ronow(i,inds(i))*Unow(i,inds(i)));
%                 ys = obj.yBL(1:inds(i));
%                 value(i) = trapz(ys, integrand);
%             end
%         end
% 
%         function value = get.Res(obj)
%             value = obj.ssurf*obj.Uinf*obj.roinf/obj.muinf;
%         end

        
    end


end