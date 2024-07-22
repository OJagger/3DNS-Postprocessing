classdef MISES_correlations
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        gam = 1.4;
    end

    methods (Static)
        function Hk = fHk(Me, H)
            Hk = (H - 0.290*Me.^2) ./ (1 + 0.113*Me.^2);
        end
    
        function H = fH(Me, Hk)
            H = 0.290*Me.^2 + Hk.*(1+0.113*Me.^2);
        end
    
        function Hs = fHs(Me, Hks)
            Hs = (Hks + 0.028*Me.^2) ./ (1 + 0.014*Me.^2);
        end
    
        function Hs = fHs2(Me, H, Rth)  % Energy shape factor straight from Me and H
            Hk = fHk(Me, H);
            Hks = fHks(Me, Rth);
            Hs = fHs(Me, Hks);
        end
        
        function Hss = fHss(Me, Hk)
            Hss = (0.064./(Hk-0.8) + 0.251).*Me.^2;
        end
        
        function Us = fUs(Hs, Hk, H)
            Us = 0.5*Hs.*(1-4*(Hk-1)./(3*H));
        end
        
        function CtEQ = fCtEQ(Hs, Us, H, Hk)
            CtEQ = Hs.*(0.015./(1-Us)).*(Hk-1).^3./(Hk.^2.*H);
        end
        
        function Ret = fRet(Me, theta, T0, P0)
            cp = 1005;
            gam = 1.4;
            R = cp*(gam-1)/gam;
            ro0 = P0/(R*T0);
            T = T0./(1+0.5*(gam-1)*Me.^2);
            ro = ro0./(1+0.5*(gam-1)*Me.^2).^(1/(gam-1));
            Ue = Me.*sqrt(gam*R*T);
            mu_ref = 5.83247e-004;
            Tref = 273;
            mu_s = 10.4;
            mue = mu_ref*(T/Tref)^(1.5)*(Tref+mu_s)/(T + mu_s);
            Ret = ro*Ue*theta/mue;
        end
        
        function Cf = fCf(Hk, Ret, Me)
            Fc = sqrt(1+0.2*Me.^2);
            Cf = (0.3*exp(-1.33*Hk).*(log10(Ret./Fc)).^(-1.74-0.31.*Hk) ...
                + 0.00011*(tanh(4-Hk/0.875) - 1))./Fc;
        end
        
        function Hks = fHks(Hk, Ret)  
            if Ret < 400
                H0 = 4;
            else
                H0 = 3+ 400./Ret;
            end
            if Hk < H0
                Hks = 1.505 + 4./Ret + (0.165-1.6./sqrt(Ret)).*(H0-Hk).^1.6./Hk;
            else
                Hks = 1.505 + 4./Ret + (Hk-H0).^2 .* (0.04./Hk + 0.007*log(Ret)./(Hk - H0 + 4./log(Ret)).^2);
            end
        end
        
    
        function Cd = fCd(Cf, Us, Ct)
            Cd = 0.5*Cf.*Us + Ct.*(1-Us);
        end
    
        function del = fDel(th, Hk, ds)
            del = th.*(3.15 + 1.72./(Hk-1)) + ds;
        end

        function Hs = fHs_ds(delStar, theta, Ma, Reth)
            H = delStar/theta;
            Hk = MISES_correlations.fHk(Ma, H);
            Hks = MISES_correlations.fHks(Hk, Reth);
            Hs = MISES_correlations.fHs(Ma, Hks);
        end

    end
end