classdef MISES_correlations
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        gam = 1.4;
    end

    

    methods (Static)
        function Hk = fHk(Me, H)
            Hk = (H - 0.290*Me.^2) ./ (1 + 0.113*Me.^2);
            if ~isreal(Hk)
                error('Hk complex')
            end
        end
    
        function H = fH(Me, Hk)
            H = 0.290*Me.^2 + Hk.*(1+0.113*Me.^2);
            if ~isreal(H)
                error('H complex')
            end
        end
    
        function Hs = fHs(Me, Hks)
            Hs = (Hks + 0.028*Me.^2) ./ (1 + 0.014*Me.^2);
            if ~isreal(Hs)
                error('Hs complex')
            end
        end
    
        function Hs = fHs2(Me, H, Rth)  % Energy shape factor straight from Me and H
            Hk = fHk(Me, H);
            Hks = fHks(Me, Rth);
            Hs = fHs(Me, Hks);
            if ~isreal(Hs)
                error('Hs complex')
            end
        end
        
        function Hss = fHss(Me, Hk)
            Hss = (0.064./(Hk-0.8) + 0.251).*Me.^2;
            if ~isreal(Hss)
                error('Hss complex')
            end
        end
        
        function Us = fUs(Hs, Hk, H)
            Us = 0.5*Hs.*(1-4*(Hk-1)./(3*H));
             % 0.5*hstar.*(1 - (4/3)*(h-1)./h);
            if ~isreal(Us)
                error('Us complex')
            end
            if any(Us < 0 | Us  > 1)
                disp('')
                % error('Us out of range');
            end
        end
        
        function CtEQ = fCtEQ(Hs, Us, H, Hk)
            CtEQ = Hs.*(0.015./(1-Us)).*(Hk-1).^3./(Hk.^2.*H);
            % (0.015*hstar./(1-us)).*( (h-1).^3 )./(h.*(h.^2))
            if ~isreal(CtEQ)
                error('CtEQ complex')
            end
        end
        
        % function Ret = fRet(Me, theta, T0, P0)
        %     cp = 1005;
        %     gam = 1.4;
        %     R = cp*(gam-1)/gam;
        %     ro0 = P0/(R*T0);
        %     T = T0./(1+0.5*(gam-1)*Me.^2);
        %     ro = ro0./(1+0.5*(gam-1)*Me.^2).^(1/(gam-1));
        %     Ue = Me.*sqrt(gam*R*T);
        %     mu_ref = 5.83247e-004;
        %     Tref = 273;
        %     mu_s = 10.4;      % Use mu_s from input file!!
        %     mue = mu_ref*(T/Tref)^(1.5)*(Tref+mu_s)/(T + mu_s);
        %     Ret = ro*Ue*theta/mue;
        %     if ~isreal(Ret)
        %         error('Ret complex')
        %     end
        % end
        
        function Cf = fCf(Hk, Ret, Me)
            Hk = reshape(Hk, [], 1);
            Ret = reshape(Ret, [], 1);
            Me = reshape(Me, [], 1);
            Fc = sqrt(1+0.2*Me.^2);
            Cf = (0.3*exp(-1.33*Hk).*(log10(Ret./Fc)).^(-1.74-0.31.*Hk) ...
                + 0.00011*(tanh(4-Hk/0.875) - 1))./Fc;
            if ~isreal(Cf)
                error('CF complex')
            end
        end

        function Hks = fHks_Hs(Hs, Me)
            Hks = (1+0.014*Me.^2).*Hs - 0.028*Me.^2;
            if ~isreal(Hks)
                error('Hks complex')
            end
        end
        
        function Hks = fHks(Hk, Ret)
            
            if isscalar(Ret)
                Ret = Ret*ones(size(Hk));
            end
            if isscalar(Hk)
                Hk = Hk*ones(size(Ret));
            end

            H0(Ret<400) = 4;
            H0(Ret >= 400) = 3+400./Ret(Ret>=400);

            Hks(Hk<H0) = 1.505 + 4./Ret(Hk<H0) + (0.165-1.6./sqrt(Ret(Hk<H0))).*(H0(Hk<H0)-Hk(Hk<H0)).^1.6./Hk(Hk<H0);
            Hks(Hk>H0) = 1.505 + 4./Ret(Hk>H0) + (Hk(Hk>H0)-H0(Hk>H0)).^2 .* (0.04./Hk(Hk>H0) + ...
                0.007*log(Ret(Hk>H0))./(Hk(Hk>H0) - H0(Hk>H0) + 4./log(Ret(Hk>H0))).^2);

            if any(~isreal(Hks))
                error('Hks complex')
            end
        end
        
        function Cd = fCd(Cf, Us, Ct)
            Cf = reshape(Cf, [], 1);
            Us = reshape(Us, [], 1);
            Ct = reshape(Ct, [], 1);
            Cd = 0.5*Cf.*Us + Ct.*(1-Us);
            if ~isreal(Cd)
                error('Cd complex')
            end
        end
        
        function del = fDel(th, Hk, ds)
            del = th.*(3.15 + 1.72./(Hk-1)) + ds;
            if ~isreal(del)
                error('del complex')
            end
        end
        
        function Hs = fHs_ds(delStar, theta, Ma, Reth)
            H = delStar/theta;
            Hk = MISES_correlations.fHk(Ma, H);
            Hks = MISES_correlations.fHks(Hk, Reth);
            Hs = MISES_correlations.fHs(Ma, Hks);
            if ~isreal(Hs)
                error('Hs complex')
            end
        end
        
        function Hk = fHk_Hks(Hks, Rt)

            for i=1:length(Hks)
                hknow = 2.0;
                Reth = Rt(i);
                H0 = 4;
                if(Reth>400)
                    H0 = 3 + (400/Reth);
                end

                del = 1e-4;

                Hksmin = 1.505 + (4)./Reth;

                if Hks(i) < Hksmin
                    hknow = H0;
                else

                    % solve for h
                    for iter=1:100
                    if(hknow<H0)   
                    hstar_guess = 1.505 + (4)./Reth + (0.165 - (1.6)./sqrt(Reth))*((H0-hknow).^1.6)./hknow;
                    else
                    hstar_guess = 1.505 + (4./Reth) + ((hknow-H0).^2).*( (0.04./hknow) + 0.007*log(Reth)./((hknow - H0 + (4./log(Reth))).^2) );    
                    end
                    if hknow < H0
                        grad = (hstar_guess - (1.505 + (4)./Reth + (0.165 - (1.6)./sqrt(Reth))*((H0-(hknow-del)).^1.6)./(hknow-del)))/del;
                        grad = sign(grad)*max(abs(grad), 0.01);
                        dh = (Hks(i) - hstar_guess)/grad;%(-0.076*( (4-hknow).^2 + 2*hknow*(4-hknow))/(hknow*hknow));
                    else
                        grad = ((1.505 + (4./Reth) + (((hknow+del)-H0).^2).*( (0.04./(hknow+del)) + 0.007*log(Reth)./(((hknow+del) - H0 + (4./log(Reth))).^2) )) - hstar_guess)/del;
                        grad = sign(grad)*max(abs(grad), 0.01);
                        dh = (Hks(i) - hstar_guess)/grad;%(0.040*( 2*hknow*(hknow-4) - (hknow-4).^2)/(hknow*hknow));
                    end
    
                    dh = dh/2^(floor(iter/10));
                    
                    if(abs(dh)<1e-4)
                        break
                    end
                    hknow = hknow + dh;
                    end
                end

                Hk(i) = hknow;
            end

            if ~isreal(Hk)
                error('Hk complex')
            end
        end

    end
end