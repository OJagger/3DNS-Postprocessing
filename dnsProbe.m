classdef dnsProbe < handle
    % DNSPROBE Class containing 3DNS probe data and metadata
    %   

    properties
        nb;
        i;
        j;
        k;
        x;
        y;
        t;
        ro;
        u;
        v;
        w;
        Et;
        vel;
        gas;
        nSkip;
        nSamples;
        T;
        Fs;
    end

    properties (Dependent = true)
        tke;
        Tu;
    end

    methods
        function obj = dnsProbe(rundir, nProbe, nSkip, blkData, gas, casetype)
            obj.nSkip = nSkip;
            obj.x = blkData.x;
            obj.y = blkData.y;
            obj.nb = blkData.nb;
            obj.i = blkData.i;
            obj.j = blkData.j;
            obj.k = blkData.k;
            obj.gas = gas;
            
            if ~isempty(rundir)
                try
                    fid = fopen(fullfile(rundir,['probe_' num2str(nProbe)]),'r');
                    switch casetype
                        case 'cpu'
                            A = fscanf(fid, '%f %f %f %f %f %f', [6 inf]);

                        case 'gpu'
                            A = fread(fid,inf,'float64');
                            A = reshape(A,6,length(A)/6);
                    end
                    fclose(fid);
                    
                    obj.t = A(1,:);
                    obj.ro = A(2,:);
                    obj.u = A(3,:)./obj.ro;
                    obj.v = A(4,:)./obj.ro;
                    obj.w = A(5,:)./obj.ro;
                    obj.Et = A(6,:);
                    
                    obj.nSamples = length(obj.t);
                catch
                    obj.t = [];
                    obj.ro = [];
                    obj.u = [];
                    obj.v = [];
                    obj.w = [];
                    obj.Et = [];

                    obj.nSamples = 0;
                end
            end
        end

        function calc_secondary_quantities(obj)

            obj.T = mean(diff(obj.t));
            obj.Fs = 1/obj.T;
%             udash = obj.u - mean(obj.u);
%             vdash = obj.v - mean(obj.v);
%             wdash = obj.w - mean(obj.w);
%             
%             obj.tke = 0.5*sqrt(mean(udash.^2)+mean(vdash.^2)+mean(wdash.^2));
%             obj.uturb = sqrt(2*obj.tke/3);
%             obj.diss = obj.uturb^3/Lturb;

        end

        function value = get.tke(obj)
            udash = obj.u-mean(obj.u);
            vdash = obj.v-mean(obj.v);
            wdash = obj.w-mean(obj.w);
            value = 0.5*(udash.^2 + vdash.^2 + wdash.^2);
        end

        function value = get.Tu(obj)
            %TURB_INTENSITY 
            %   Calculate turbulence intensuty
            
            U = mean(sqrt(obj.u.^2 + obj.v.^2 + obj.w.^2));
            value = sqrt(2*mean(obj.tke)/3)/U;

        end

        function value = plot_spectrum(obj,ax,plot_inertial, fmt)

            if nargin < 4 || isempty(fmt)
                fmt = [];
            end 
            if nargin < 3 || isempty(plot_inertial)
                plot_inertial = false;
            end
            
            if nargin<2 || isempty(ax)
                ax = gca;
            end


            obj.calc_secondary_quantities;

            udash = obj.u-mean(obj.u);
            vdash = obj.v-mean(obj.v);
            wdash = obj.w-mean(obj.w);

            numlags = floor(length(udash)/6);

            [R_u, lags_u] = autocorr(obj.u,'NumLags',numlags);
            [R_v, lags_v] = autocorr(obj.v,'NumLags',numlags);
            [R_w, lags_w] = autocorr(obj.w,'NumLags',numlags);

            T_int_u = trapz(obj.t(1:numlags+1), R_u);
            T_int_v = trapz(obj.t(1:numlags+1), R_v);
            T_int_w = trapz(obj.t(1:numlags+1), R_w);

            T_int = max([T_int_u T_int_v T_int_w]);
            windowL = 500 % floor(T_int * obj.Fs);
            [PSD, F_psd] = pwelch(sqrt(udash.^2 + vdash.^2 + wdash.^2), ...
                windowL, floor(windowL/2), [], obj.Fs, 'psd');

            

%             N = length(obj.tke);
%             Ts = (obj.t(end)-obj.t(1))/N;
%             fs = 1/Ts;
%             Y = fft(obj.tke);
%             P2 = abs(Y/N);
%             P1 = P2(1:N/2+1);
%             P1(2:end-1) = 2*P1(2:end-1);
%             psd = P1/fs;
%             f = fs*(0:(N/2))/N;
% 
            value = loglog(F_psd, PSD, fmt);
            if plot_inertial
                B = 50000; % B is a fudge to adjust the y intercept
                fit_wp = logspace(log10(min(F_psd(2:end))),6,100);
                fit_Epw = B.*fit_wp.^(-5/3);
                hold on
                loglog(fit_wp, fit_Epw, 'k:')
            end
            
            grid on
            xlabel('f (Hz)')
            ylabel('PSD (dB/Hz)')
            set(gca,'fontSize',12)
        end

        function concatenate(obj, newProbe)
            if obj.i ~= newProbe.i; fprintf('i mismatch\n'); return; end
            if obj.j ~= newProbe.j; fprintf('i mismatch\n'); return; end
            if obj.nb ~= newProbe.nb; fprintf('block mismatch\n'); return; end
            if obj.nSkip ~= newProbe.nSkip; fprintf('nSkip mismatch\n'); return; end

            obj.u = [obj.u newProbe.u];
            obj.v = [obj.v newProbe.v];
            obj.w = [obj.w newProbe.w];
            obj.ro = [obj.ro newProbe.ro];
            obj.Et = [obj.Et newProbe.Et];
            obj.t = [obj.t newProbe.t];
            obj.nSamples = obj.nSamples + newProbe.nSamples;

        end

    end
end