function ydash = MISES_Newton_BL(xsurf, Mesurf, Uesurf, gas, bcs, th0, delStar0, Ctau0)

    if nargin < 7 || isempty(delStar0)
        H = 1.4;
        delStar0 = H*th0;
    end

    xsurf = xsurf-xsurf(1);

    h = 1e-5*[th0(1); th0(1); delStar0(1); delStar0(1); Ctau0(1); Ctau0(1)];
%     if nargin < 8 || isempty(Ctau0)
%         Me = 1.2;
%         M
%         Ctau0 = 


%     global x1 x2 Ue1 Ue2 M1 M2 nua

    N = length(xsurf);
    Q(1:N, 1) = th0;
    Q(1:N, 2) = delStar0;
    Q(1:N, 3) = Ctau0;

    rgas = gas.cp*(gas.gam-1)/gas.gam;
    T = bcs.Toin*T_T0(Mesurf, gas.gam);
    rosurf = (bcs.Poin/rgas/bcs.Toin)*ro_ro0(Mesurf, gas.gam);
    nusurf = sutherland_mu(T, gas.mu_ref, gas.mu_cref, gas.mu_tref)./rosurf;
    
%     Me = interp1(xsurf, Mesurf, x);
%     Uesurf = Mesurf * interp1(xsurf, Uesurf, x);

    % Residual vector
    F = zeros(3*N,1);

    % Jacobian matrix;
    J = zeros(3*N, 3*N);

    % Jacobian functor
    % Jac = @(x,h,fun)((fun(repmat(x,size(x'))+diag(h))-fun(repmat(x,size(x'))))./(100*h)');
    th_fun = @(u) (R_theta(u));
    H_fun = @(u) (R_Hs(u));
    ct_fun = @(u) (R_Ct(u));

    % Your function
%     f = @(x)[x(1)^2 + x(2)^2; x(1)^3.*x(2)^3];
%     % Point at which to estimate it
%     x = [1;1];
%     % Step to take on each dimension (has to be small enough for precision)
%     h = 1e-5*ones(size(x));
%     % Compute the jacobian
%     J(x,h,f)

    for i=1:N-1
        x1 = xsurf(i);
        x2 = xsurf(i+1);
        Ue1 = Uesurf(i);
        Ue2 = Uesurf(i+1);
        M1 = Mesurf(i);
        M2 = Mesurf(i+1);
        nua = 0.5*(nusurf(i)+nusurf(i+1));

        xv = reshape(Q(i:i+1, :), [], 1);
        F(3*(i-1)+1) = R_theta(xv);
        F(3*(i-1)+2) = R_Hs(xv);
        F(3*(i-1)+3) = R_Ct(xv);

        J(3*(i-1)+1, 3*(i-1)+1:3*(i-1)+6) = Jac(xv, h, th_fun);
        J(3*(i-1)+2, 3*(i-1)+1:3*(i-1)+6) = Jac(xv, h, H_fun);
        J(3*(i-1)+3, 3*(i-1)+1:3*(i-1)+6) = Jac(xv, h, ct_fun);

    end

    % J = real(J);
    F(end-2) = th0(1) - Q(1,1);
    F(end-1) = delStar0(1) - Q(1,2);
    F(end) = Ctau0(1) - Q(1,3);

    cond(J)
    a = 1;
    J(end-2, 1) = a;
    J(end-1, 2) = a;
    J(end, 3) = a;
    cond(J)
    

    dQ = J\(-F);
    dQ = reshape(dQ, 3, N)';
    %dQ = max()
    Q = Q+dQ;

    function J = Jac(x, h, fun)
        xm = repmat(x,size(x'));
        f0 = fun(xm);
        f1 = fun(xm+diag(h));
        J = ((f1-f0))./(h)';
    end
    

    function R = R_theta(xv) % Calculate residual of momentum integral equation

        nj = size(xv, 2);
        R = zeros(1, nj);

        for j=1:nj
            th1 = xv(1,j);
            th2 = xv(2,j);
            ds1 = xv(3,j);
            ds2 = xv(4,j);
            ct1 = xv(5,j);
            ct2 = xv(6,j);
            
            xa = 0.5*(x1+x2);
            tha = 0.5*(th1+th2);
            dsa = 0.5*(ds1+ds2);
            cta = 0.5*(ct1+ct2);
            Uea = 0.5*(Ue1+Ue2);
            Ma = 0.5*(M1+M2);
            
            Reta = tha*Uea/nua;
            Ha = dsa/tha;
            Hka = fHk(Ma, Ha);
            Cf = fCf(Hka, Reta, Ma);
    
            R(j) = log(th2/th1)/log(x2/x1) - (xa/tha)*Cf/2 + (Ha + 2 - Ma^2) * log(Ue2/Ue1)/log(x2/x1);
%             R(j) = (th2-th1)/(x2-x1) -Cf/2 + (Ha + 2 - Ma^2) * (Ue2-Ue1)*tha/Uea/(x2-x1);
        end
        R;
    end

    function R = R_Hs(xv) % Calculate residual of shape factor equation

        nj = size(xv, 2);
        R = zeros(1, nj);

        for j=1:nj
            th1 = xv(1,j);
            th2 = xv(2,j);
            ds1 = xv(3,j);
            ds2 = xv(4,j);
            ct1 = xv(5,j);
            ct2 = xv(6,j);
    
            xa = 0.5*(x1+x2);
            tha = 0.5*(th1+th2);
            dsa = 0.5*(ds1+ds2);
            cta = 0.5*(ct1+ct2);
            Uea = 0.5*(Ue1+Ue2);
            xa = 0.5*(x1+x2);
            Ma = 0.5*(M1+M2);
    
            Reta = tha*Uea/nua;
            Ret1 = th1*Ue1/nua;
            Ret2 = th2*Ue2/nua;
            Ha = dsa/tha;
            Hka = fHk(Ma, Ha);
            Hks = fHks(Ma, Reta);
            Hs = fHs(Ma, Hks);
            Cf = fCf(Hka, Reta, Ma);
            Hs1 = fHs2(M1, ds1/th1, Ret1);
            Hs2 = fHs2(M2, ds2/th2, Ret2);
            Us = fUs(Hs, Hka, Ha);
            Cd = fCd(Cf, Us, cta*cta);
            Hss = fHss(Ma, Hka);
    
            R(j) = log(Hs2/Hs1)/log(x2/x1) + (xa/tha)*(Cf/2 - 2*Cd/Hs) + (2*Hss/Hs + 1 - Ha)*log(Ue2/Ue1)/log(x2/x1);        
%             R(j) = (Hs2-Hs1)/(x2-x1) + (Hs*Cf/2 - 2*Cd)/tha + (2*Hss + Hs*(1 - Ha))*(Ue2-Ue1)/(x2-x1)/Uea;     
        end
    end

    
    function R = R_Ct(xv)       % Residual of shear lag equation

        nj = size(xv, 2);
        R = zeros(1, nj);

        for j=1:nj
            th1 = xv(1,j);
            th2 = xv(2,j);
            ds1 = xv(3,j);
            ds2 = xv(4,j);
            ct1 = xv(5,j);
            ct2 = xv(6,j);
    
            tha = 0.5*(th1+th2);
            dsa = 0.5*(ds1+ds2);
            Ma = 0.5*(M1+M2);
            Uea = 0.5*(Ue1+Ue2);
    
            Reta = tha*Uea/nua;
            Rth2 = th2*Ue2/nua;
            Ha = dsa/tha;
            H2 = ds2/th2;
            Hka = fHk(Ma, Ha);
            Hk2 = fHk(M2, H2);
            Hks = fHks(Ma, Reta);
            Hs = fHs(Ma, Hks);
            Us = fUs(Hs, Hka, Ha);
            CtEQ = sqrt(fCtEQ(Hs, Us, Ha, Hka));
            del = fDel(th2, Hk2, ds2);
    
            R(j) = 2*(del/ct2)*(ct2-ct1)/(x2-x1) - 4.2*(CtEQ - ct2);
        end
    end

    function Hk = fHk(Me, H)
        Hk = (H - 0.290*Me^2) / (1 + 0.113*Me^2);
    end

    function H = fH(Me, Hk)
        H = 0.290*Me^2 + Hk*(1+0.113*Me^2);
    end

    function Hs = fHs(Me, Hks)
        Hs = (Hks + 0.028*Me^2) / (1 + 0.014*Me^2);
    end

    function Hs = fHs2(Me, H, Rth)  % Energy shape factor straight from Me and H
        Hk = fHk(Me, H);
        Hks = fHks(Me, Rth);
        Hs = fHs(Me, Hks);
    end
    
    function Hss = fHss(Me, Hk)
        Hss = (0.064/(Hk-0.8) + 0.251)*Me^2;
    end
    
    function Us = fUs(Hs, Hk, H)
        Us = 0.5*Hs*(1-4*(Hk-1)/(3*H));
    end
    
    function CtEQ = fCtEQ(Hs, Us, H, Hk)
        CtEQ = Hs*(0.015/(1-Us))*(Hk-1)^3/(Hk^2*H);
    end
    
    function Ret = fRet(Me, theta, T0, P0)
        cp = 1005;
        gam = 1.4;
        R = cp*(gam-1)/gam;
        ro0 = P0/(R*T0);
        T = T0/(1+0.5*(gam-1)*Me^2);
        ro = ro0/(1+0.5*(gam-1)*Me^2)^(1/(gam-1));
        Ue = Me*sqrt(gam*R*T);
        mu_ref = 5.83247e-004;
        Tref = 273;
        mu_s = 10.4;
        mue = mu_ref*(T/Tref)^(1.5)*(Tref+mu_s)/(T + mu_s);
        Ret = ro*Ue*theta/mue;
    end
    
    function Cf = fCf(Hk, Ret, Me)
        Fc = sqrt(1+0.2*Me^2);
        Cf = (0.3*exp(-1.33*Hk)*(log10(Ret/Fc))^(-1.74-0.31*Hk) ...
            + 0.00011*(tanh(4-Hk/0.875) - 1))/Fc;
    end
    
    function Hks = fHks(Hk, Ret)  
        if Ret < 400
            H0 = 4;
        else
            H0 = 3+ 400/Ret;
        end
        if Hk < H0
            Hks = 1.505 + 4/Ret + (0.165-1.6/sqrt(Ret))*(H0-Hk)^1.6/Hk;
        else
            Hks = 1.505 + 4/Ret + (Hk-H0)^2 * (0.04/Hk + 0.007*log(Rth)/(Hk - H0 + 4/log(Rth))^2);
        end
    end

    function Cd = fCd(Cf, Us, Ct)
        Cd = 0.5*Cf*Us + Ct*(1-Us);
    end

    function del = fDel(th, Hk, ds)
        del = th*(3.15* + 1.72/(Hk-1)) + ds;
    end

end