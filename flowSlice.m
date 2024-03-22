classdef flowSlice < handle
    %FLOWSLICE Generic class containing a slice of the flow

    properties
        NB;
        gas;
        bcs;
        ro;
        u;
        v;
        w;
        Et;
%         xSurf;
%         yBL;
%         xO;
%         yO;
%         iO;
%         jO;
%         blkO;
%         oblocks;
%         oblocks_flip;
        blk;
%         iLE;
%         iTE;
%         n;
%         ssurf;          % Surface distance fron LE
%         vortZ;          % Z vorticity
    end

    properties (Dependent = true, Hidden = true)
        T;             % Temperature
        p;             % p stat
        StR;            % Strain rate magnitude
        M;              % Mach No
        s;              % Entropy ( cp*log(T/300) - R*log(p/1e5) )
        sfunc;
        ros;            % Entropy per unit volume
        vel;            % Velocity
        mu;             % Viscosity
        nu;             % Kinematic viscosity
        p0;
        T0;
        schlieren;      % |grad(ro)|/ro
        cellSize;
        St;
        St_an              % Traceless strain
        S_an_mag;       % Magnitude of anisotropic componant of strain tensor
        local_cfl;
        wallDist;       % Distance from wall
        flowAng;
        muSij2;         % Dissipation due to mean strain
    end

    methods
        function obj = flowSlice(blk, gas, bcs)
            %FLOWSLICE Construct a flowSlice object
            disp('Constructing flowSlice')

            obj.gas = gas;
            obj.bcs = bcs;
            obj.NB = size(blk.blockdims,1);
            obj.blk = blk;
        end         % End of constructor

        function value = get.p(obj)
            value = obj.get_p;
        end

        function obj = set.p(obj, value)
            obj.set_p(value);
        end

        function value = get_p(obj)
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                value{nb} = (obj.gas.gam -1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
            end
        end

        function value = get.T(obj)
            value = obj.get_T;
        end

        function obj = set.T(obj, value)
            obj.set_T(value);
        end

        function value = get.ros(obj)
            value = obj.get_ros;
        end

        function obj = set.ros(obj, value)
            obj.set_ros(value)
        end

        function value = get.StR(obj)
            value = obj.get_StR;
        end
        

        function obj = set.StR(obj, value)
            obj.set_StR(value);
        end

        function set_p(obj,value)
        end

        function set_T(obj,value)
        end

        function set_StR(obj,value)
        end

        function set_ros(obj,value)
        end

        function value = get_T(obj)
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                value{nb} = obj.p{nb}./(obj.ro{nb}*obj.gas.rgas);
            end
        end


        function value = get_ros(obj)
            value = cell(1,obj.NB);
            Tnow = obj.T;
            pnow = obj.p;
            for nb = 1:obj.NB
                value{nb} = obj.ro{nb}.*(obj.gas.cp*log(Tnow{nb}/obj.bcs.Toin) - ...
                    obj.gas.cp*(1-1/obj.gas.gam)*log(pnow{nb}/obj.bcs.Poin));
            end
        end

        function value = get.vel(obj)
            disp('Calculating vel')
            value = cell(1,obj.NB);
            for nb = 1:obj.NB
                value{nb} = sqrt(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2);
            end
        end

        function value = get.p0(obj)
            pnow = obj.p;
            Mnow = obj.M;
            for nb = 1:obj.NB
                value{nb} = pnow{nb}.*(1+0.5*(obj.gas.gam-1)*Mnow{nb}.^2).^(obj.gas.gam/(obj.gas.gam-1));
            end
        end

        function value = get.T0(obj)
            tnow = obj.T;
            Mnow = obj.M;
            for nb = 1:obj.NB
                value{nb} = tnow{nb}.*(1+0.5*(obj.gas.gam-1)*Mnow{nb}.^2);
            end
        end

        function value = get.M(obj)
            disp('Calculating M')
            value = cell(1,obj.NB);
            pnow = obj.p;
            for nb = 1:obj.NB
                %pnow = (obj.gas.gam - 1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                Tnow = pnow{nb}./(obj.ro{nb}*obj.gas.rgas);
                velnow = sqrt(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2);
                value{nb} = velnow./sqrt(obj.gas.gam*obj.gas.rgas*Tnow);
            end
        end

        function value = get.s(obj)
            disp('Calculating s')
            value = cell(1,obj.NB);
            for nb = 1:obj.NB
                pnow = (obj.gas.gam - 1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                Tnow = pnow./(obj.ro{nb}*obj.gas.rgas);
                value{nb} = obj.gas.cp*log(Tnow/obj.bcs.Toin) - obj.gas.rgas*log(pnow/obj.bcs.Poin);
            end
        end

        function value = get.sfunc(obj)
            value = cell(1,obj.NB);
            for nb = 1:obj.NB
                pnow = (obj.gas.gam - 1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                Tnow = pnow./(obj.ro{nb}*obj.gas.rgas);
                value{nb} = obj.gas.cp*log(Tnow/obj.bcs.Toin) - obj.gas.rgas*log(pnow/obj.bcs.Poin);
            end
        end

        function value = get.mu(obj)
            disp('Calcualting mu')
            value = cell(1,obj.NB);
            for nb = 1:obj.NB
                pnow = obj.p{nb};%(obj.gas.gam - 1)*(obj.Et{nb} - 0.5*(obj.u{nb}.^2 + obj.v{nb}.^2 + obj.w{nb}.^2).*obj.ro{nb});
                Tnow = pnow./(obj.ro{nb}*obj.gas.rgas);
                value{nb} = sutherland_mu(Tnow, obj.gas.mu_ref, obj.gas.mu_cref, obj.gas.mu_tref);%obj.gas.mu_ref*(Tnow/obj.gas.mu_tref).^(3/2) .* (obj.gas.mu_cref + obj.gas.mu_tref)./(obj.gas.mu_cref + Tnow);
            end
        end

        function value = get.flowAng(obj)
            disp('Calculating flow angle')
            value = cell(1,obj.NB);
            for nb = 1:obj.NB
                value{nb} = atan2d(obj.v{ib}./obj.u{ib});
            end
        end

        function value = get.nu(obj)
            munow = obj.mu;
            ronow = obj.ro;
            for nb = 1:obj.NB
                value{nb} = munow{nb}./ronow{nb};
            end
        end

        

        function value = get.cellSize(obj)
            fprintf('Calculating Cell Sizes\n')
            dz = obj.blk.span/(obj.blk.nk{1}-1);
            
            value = {};
            for ib = 1:obj.NB
                ni = size(obj.blk.x{ib},1);
                nj = size(obj.blk.x{ib},2);
                area = zeros(ni-1, nj-1);
                for i=1:ni-1
                    for j=1:nj-1
                        xnow = [obj.blk.x{ib}(i,j) obj.blk.x{ib}(i+1,j) ...
                            obj.blk.x{ib}(i+1,j+1) obj.blk.x{ib}(i,j+1)];
                        ynow = [obj.blk.y{ib}(i,j) obj.blk.y{ib}(i+1,j) ...
                            obj.blk.y{ib}(i+1,j+1) obj.blk.y{ib}(i,j+1)];
                        area(i,j) = abs(polyarea(xnow,ynow));
                    end
                end
                area = dz*area;
                value{ib}(1,1) = area(1,1);
                value{ib}(1,nj) = area(1,nj-1);
                value{ib}(ni,1) = area(ni-1,1);
                value{ib}(ni,nj) = area(ni-1,nj-1);
                for i = 2:ni-1
                    value{ib}(i,1) = 0.5*(area(i-1,1)+area(i,1));
                    value{ib}(i,end) = 0.5*(area(i-1,end)+area(i,end));
                end
                for j = 2:nj-1
                    value{ib}(1,j) = 0.5*(area(1,j-1)+area(1,j));
                    value{ib}(end,j) = 0.5*(area(end,j-1)+area(end,j));
                end
                for i=2:ni-1
                    for j=2:nj-1
                        value{ib}(i,j) = 0.25*(area(i-1,j-1)+area(i-1,j)+area(i,j-1) +area(i,j));
                    end
                end
                value{ib} = value{ib}.^(1/3);
            end
        end


        function value = get.schlieren(obj)
            disp('calculating grad(ro)/ro')
            value = cell(1,obj.NB);
            for nb=1:obj.NB
                [drodx, drody] = gradHO(obj.blk.x{nb},obj.blk.y{nb},obj.ro{nb});
                value{nb} = sqrt(drodx.^2 + drody.^2)./obj.ro{nb};
            end
        end

        function getSize(obj)
            props = properties(obj);
            totSize = 0; 
            for ii=1:length(props) 
                currentProperty = obj.(props{ii});
                temp = whos('currentProperty'); 
                totSize = totSize + temp.bytes; 
            end
          
            fprintf(1, '%d MB\n', totSize/1e6);
        end

        function value = get_StR(obj)
            value = cell(1,obj.NB);
            for nb =1:obj.NB
                value{nb} = strain_rate_magnitude(obj.blk.x{nb}, obj.blk.y{nb}, obj.u{nb}, obj.v{nb});
            end
        end
        

        function value = get.St_an(obj)

            S = obj.St;

            for ib = 1:obj.NB

                tr = S{ib}(:,:,1,1) + S{ib}(:,:,2,2) + S{ib}(:,:,3,3);

                St = S{ib};
                for i=1:3
                    St(:,:,i,i) = St(:,:,i,i) - tr/3;
                end

                value{ib} = St;

%                 [DUDX,DUDY] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.u{ib});
%                 [DVDX,DVDY] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.v{ib});
% 
%                 %Traceless strain tensor
%                 S = zeros(obj.blk.blockdims(ib,1),obj.blk.blockdims(ib,2),3,3);
%                 
%                 S(:,:,1,1) = 2*DUDX/3 - DVDY/3;
%                 S(:,:,2,2) = 2*DVDY/3 - DUDX/3;
%                 S(:,:,3,3) = -(DUDX+DVDY)/3;
% 
%                 S(:,:,1,2) = 0.5*(DUDY+DVDX);
%                 S(:,:,2,1) = S(:,:,1,2);

            end
        end

        function value = St_an2(obj)

            S = obj.St_an;
            for ib = 1:obj.NB
                value{ib} = sum(sum(S{ib}.*S{ib},4),3);
            end
        end

        function value = get.St(obj)

            for ib = 1:obj.NB

                [DUDX,DUDY] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.u{ib});
                [DVDX,DVDY] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.v{ib});

                %Traceless strain tensor
                S = zeros(obj.blk.blockdims(ib,1),obj.blk.blockdims(ib,2),3,3);
                
                S(:,:,1,1) = DUDX;
                S(:,:,2,2) = DVDY;

                S(:,:,1,2) = 0.5*(DUDY+DVDX);
                S(:,:,2,1) = S(:,:,1,2);

                value{ib} = S;
            end

        end


        function value = get.S_an_mag(obj)
            Snow = obj.St_an;
            for ib = 1:obj.NB
                value{ib} = sqrt(sum(sum(Snow{ib}.*Snow{ib},4),3));
            end
        end

        function value = get.muSij2(obj)

            for ib = 1:obj.NB

                [DUDX,DUDY] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.u{ib});
                [DVDX,DVDY] = gradHO(obj.blk.x{ib},obj.blk.y{ib},obj.v{ib});

                % Strain tensor
                S = zeros(obj.blk.blockdims(ib,1),obj.blk.blockdims(ib,2),3,3);
                
                S(:,:,1,1) = DUDX;
                S(:,:,2,2) = DVDY;

                S(:,:,1,2) = 0.5*(DUDY+DVDX);
                S(:,:,2,1) = S(:,:,1,2);

                value{ib} = obj.mu{ib}.*sum(sum(S.*S,4),3);
            end


        end

        function value = get.local_cfl(obj)
            cfl = 1.0;
            Vnow = obj.vel;
            Tnow = obj.T;
            c = {};
            cmax = [];
            for ib = 1:obj.NB
                c{ib} = sqrt(obj.gas.gam*obj.gas.rgas*Tnow{ib});
                cmax(ib) = max(c{ib},[],'all');
                vmax(ib) = max(Vnow{ib},[],'all');
                [dxi, dxj] = gradHOij(obj.blk.x{ib});
                [dyi, dyj] = gradHOij(obj.blk.y{ib});
                dcelli = sqrt(dxi.^2 + dyi.^2);
                dcellj = sqrt(dxj.^2 + dyj.^2);
                dcell{ib} = min(dcelli, dcellj);
            end
            vmax = max(vmax);
            cmax = max(cmax);

            for ib = 1:obj.NB
                dt(ib) = min(cfl*dcell{ib}/(vmax+cmax),[],'all');
            end
            dt = min(dt);
            fprintf('Slice Δt = %5.3e\n',dt)
            for ib = 1:obj.NB
                value{ib} = dt*(Vnow{ib}+c{ib})./dcell{ib};
            end 
        end

        function value = get.wallDist(obj)
            value = obj.blk.walldist;
        end

        function kPlot(obj, prop)
            q = obj.(prop);
            hold on
            ax = gca;
            for i=1:obj.NB
                s = pcolor(ax, obj.blk.x{i}, obj.blk.y{i}, q{i});
            end
            shading('interp');
            axis equal
        end

    end
end
